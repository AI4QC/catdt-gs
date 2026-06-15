from ase import units
from ase.io import Trajectory, read, write
from ase import Atoms, units
from ase.calculators.calculator import Calculator
from ase.constraints import FixAtoms
from ase import units
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution,Stationary,ZeroRotation
from ase.md import MDLogger
from ase.calculators.plumed import Plumed
from ase import units

from mace.calculators import MACECalculator

import integrator as md_integrator

import gc
import random
import torch
import time
import os
import yaml
import copy
from tqdm import tqdm
from pathlib import Path
import numpy as np

def merge_dicts(dict1: dict, dict2: dict):
    if not isinstance(dict1, dict):
        raise ValueError(f"Expecting dict1 to be dict, found {type(dict1)}.")
    if not isinstance(dict2, dict):
        raise ValueError(f"Expecting dict2 to be dict, found {type(dict2)}.")

    return_dict = copy.deepcopy(dict1)
    duplicates = []

    for k, v in dict2.items():
        if k not in dict1:
            return_dict[k] = v
        else:
            if isinstance(v, dict) and isinstance(dict1[k], dict):
                return_dict[k], duplicates_k = merge_dicts(dict1[k], dict2[k])
                duplicates += [f"{k}.{dup}" for dup in duplicates_k]
            else:
                return_dict[k] = dict2[k]
                duplicates.append(k)

    return return_dict, duplicates

def load_config(path: str, previous_includes: list = []):
    if path in previous_includes:
        raise ValueError(
            f"Cyclic config include detected. {path} included in sequence {previous_includes}."
        )
    previous_includes = previous_includes + [path]

    direct_config = yaml.safe_load(open(path, "r"))

    # Load config from included files.
    if "includes" in direct_config:
        includes = direct_config.pop("includes")
    else:
        includes = []
    if not isinstance(includes, list):
        raise AttributeError(
            "Includes must be a list, '{}' provided".format(type(includes))
        )

    config = {}
    duplicates_warning = []
    duplicates_error = []

    for include in includes:
        include_config, inc_dup_warning, inc_dup_error = load_config(
            include, previous_includes
        )
        duplicates_warning += inc_dup_warning
        duplicates_error += inc_dup_error

        # Duplicates between includes causes an error
        config, merge_dup_error = merge_dicts(config, include_config)
        duplicates_error += merge_dup_error

    # Duplicates between included and main file causes warnings
    config, merge_dup_warning = merge_dicts(config, direct_config)
    duplicates_warning += merge_dup_warning

    return config, duplicates_warning, duplicates_error


class AverageForceCalculator(Calculator):
    implemented_properties = ['forces', 'energy', 'potential']

    def __init__(self, calculators, **kwargs):
        super().__init__(**kwargs)
        self.calculators = calculators

    def calculate(self, atoms=None, properties=['forces', 'energy', 'potential'], system_changes=None):
        super().calculate(atoms, properties, system_changes)

        total_forces = 0
        total_energy = 0
        all_forces = []
        total_mu = 0
        all_mu = []

        for calc in self.calculators:
            atoms_data = copy.deepcopy(atoms)
            atoms_data.set_calculator(calc)
            calc.calculate(atoms_data, properties, system_changes)
            total_forces += calc.results['forces']
            total_energy += calc.results['energy']
            total_mu += calc.results['potential']
            all_forces.append(calc.results['forces'])
            all_mu.append(calc.results['potential'])


        average_forces = total_forces / len(self.calculators)
        average_energy = total_energy / len(self.calculators)
        average_mu = total_mu / len(self.calculators)
        all_forces_array = np.array(all_forces) # shape: (num_calculators, num_atoms, 3)

        # Calculate the standard deviation of the forces in each direction
        force_std = np.std(all_forces_array, axis=0) # shape: (num_atoms, 3)
        mu_std = np.std(np.array(all_mu))

        # Calculate the total standard deviation for each atom
        total_std = np.sqrt(np.sum(force_std**2, axis=1)) # shape: (num_atoms,)

        # Find the maximum total standard deviation across all atoms
        max_total_std = np.max(total_std)
        max_std_atom_index = np.argmax(total_std)

        self.results['forces'] = average_forces
        self.results['energy'] = average_energy
        self.results['mu'] = average_mu
        atoms.info['potential'] = average_mu
        self.results['force_std'] = max_total_std
        self.results['mu_std'] = mu_std
        print('force_std: ', max_total_std)
        print('fermi_std:', mu_std)

    def get_max_std(self):
        return self.results.get('force_std', 0)

    def get_mu_std(self):
        return self.results.get('mu_std', 0)

    def get_mu(self, atoms=None):
        return self.results.get('mu', 0)

class NeuralMDLogger(MDLogger):
    def __init__(self,
                 *args,
                 start_time=0,
                 verbose=True,
                 **kwargs):
        if start_time == 0:
            header = True
        else:
            header = False
        super().__init__(header=header, *args, **kwargs)
        """
        Logger uses ps units.
        """
        self.start_time = start_time
        self.verbose = verbose
        if verbose:
            print(self.hdr)
        self.natoms = self.atoms.get_number_of_atoms()

    def __call__(self):
        if self.start_time > 0 and self.dyn.get_time() == 0:
            return
        epot = self.atoms.get_potential_energy()
        ekin = self.atoms.get_kinetic_energy()
        temp = ekin / (1.5 * units.kB * self.natoms)
        if self.peratom:
            epot /= self.natoms
            ekin /= self.natoms
        if self.dyn is not None:
            t = self.dyn.get_time() / (1000*units.fs) + self.start_time
            dat = (t,)
        else:
            dat = ()
        dat += (epot+ekin, epot, ekin, temp)
        if self.stress:
            dat += tuple(self.atoms.get_stress() / units.GPa)
        self.logfile.write(self.fmt % dat)
        self.logfile.flush()

        if self.verbose:
            print(self.fmt % dat)


class Simulator:
    def __init__(self,
                 atoms,
                 integrator,
                 T_init,
                 start_time=0,
                 save_dir='./log',
                 restart=False,
                 save_frequency=100,
                 min_temp=0.1,
                 max_temp=100000):
        self.atoms = atoms
        self.integrator = integrator
        self.save_dir = Path(save_dir)
        self.min_temp = min_temp
        self.max_temp = max_temp
        self.natoms = self.atoms.get_number_of_atoms()

        # intialize system momentum
        if not restart:
            MaxwellBoltzmannDistribution(self.atoms, temperature_K = T_init, rng=np.random.RandomState(11300402))
            Stationary(self.atoms), ZeroRotation(self.atoms)

        # attach trajectory dump
        self.traj = Trajectory(self.save_dir / 'atoms.traj', 'a', self.atoms)
        self.integrator.attach(self.traj.write, interval=save_frequency)

        # attach log file
        self.integrator.attach(NeuralMDLogger(self.integrator, self.atoms,
                                        self.save_dir / 'thermo.log',
                                        start_time=start_time, mode='a'),
                               interval=save_frequency)
    def run(self, steps, max_std_threshold, max_mu_threshold):

        for step in tqdm(range(steps)):

            self.integrator.run(1)
            if step % 100 == 0:
               gc.collect()
               torch.cuda.empty_cache()

            max_std = self.atoms.get_calculator().calc.get_max_std()
            mu_std = self.atoms.get_calculator().calc.get_mu_std()

            if 0.25 > max_std > max_std_threshold:
                print(f"max_std {max_std} is greater than threshold {max_std_threshold}. Collecting structure.")
                self.atoms.write(self.save_dir / 'structure_force.xyz',append=True)

            elif mu_std > max_mu_threshold:

                print(f"mu_std {mu_std} is greater than threshold {max_mu_threshold}. Collecting structure.")
                self.atoms.write(self.save_dir / 'structure_fermi.xyz',append=True)

            ekin = self.atoms.get_kinetic_energy()
            temp = ekin / (1.5 * units.kB * self.natoms)
            if temp < self.min_temp or temp > self.max_temp:
                print(f'Temperature {temp:.2f} is out of range: \
                        [{self.min_temp:.2f}, {self.max_temp:.2f}]. \
                        Early stopping the simulation.')
                break

        self.traj.close()
        return True, step + 1


if __name__ == '__main__':

#   torch.use_deterministic_algorithms(True)

    calculator1 = MACECalculator(model_paths=['AuORR_run-1130.model'], device='cuda')
    calculator2 = MACECalculator(model_paths=['AuORR_run-2008.model'], device='cuda')
    #calculator3 = MACECalculator(model_paths=['AuORR_alkal_3_compiled.model'], device='cuda')
    #calculator4 = MACECalculator(model_paths=['AuORR_alkal_4_compiled.model'], device='cuda')
    average_calculator = AverageForceCalculator([calculator1, calculator2])

    #print(average_calculator)

    torch.set_default_dtype(torch.float64)
    config, _, _ = load_config("inputs.yml")
    #print(config)

    atoms = read("init.xyz")
    atoms.pbc = [True, True, True]
    simulated_time = 0
    simulated_step = 0
    
    atoms.set_calculator(average_calculator)
    plumed_inp = ["UNITS LENGTH=A TIME=fs ENERGY=eV",
                  "dOO: DISTANCE ATOMS=93,94",
                  "surfAu: COM ATOMS=188-203",
                  "dAuOH: DISTANCE ATOMS=93,surfAu COMPONENTS",
                  "allH: GROUP ATOMS=1-92",
                  "etotal: ENERGY",
                  "Nproton: COORDINATION GROUPA=127 GROUPB=allH R_0=1.20 D_0=0 NN=9 MM=18 NLIST NL_CUTOFF=10 NL_STRIDE=5",
                  "ghosta: FIXEDATOM AT=0,0,20",
                  "H_proton: DISTANCE ATOMS=ghosta,127 COMPONENTS",
                  "abs_H: MATHEVAL ARG=H_proton.z FUNC=abs(x) PERIODIC=NO",
                  "abs_dAuOH: MATHEVAL ARG=dAuOH.z FUNC=abs(x) PERIODIC=NO",
                  "CNendOH: COORDINATION GROUPA=93 GROUPB=allH R_0=1.20 D_0=0 NN=9 MM=18 NLIST NL_CUTOFF=10 NL_STRIDE=5",
                  "cOadH: COORDINATION GROUPA=94 GROUPB=allH R_0=1.20 D_0=0 NN=9 MM=18 NLIST NL_CUTOFF=10 NL_STRIDE=5",
                  "LOWER_WALLS LABEL=lwall ARG=Nproton,dOO AT=2.70,1.35 KAPPA=1,1",
                  "UPPER_WALLS LABEL=uwall ARG=abs_H,dOO AT=2.50,3.2 KAPPA=1,10",
		  "METAD LABEL=metad ARG=dOO SIGMA=0.20 HEIGHT=0.2 PACE=250 BIASFACTOR=15 GRID_MIN=1.1 GRID_MAX=3.6 GRID_BIN=500 CALC_RCT FILE=HILLS",
                  "PRINT ARG=dOO,abs_dAuOH,cOadH,CNendOH,Nproton,abs_H,etotal,metad.bias,metad.rbias,metad.rct,metad.work,uwall.bias,lwall.bias STRIDE=1 FILE=COLVAR",
                  "FLUSH STRIDE=1"]

    atoms.calc = Plumed(calc=average_calculator, input=plumed_inp,
                        timestep=1.0*units.fs, atoms=atoms, kT=0.02585199101165164)
    os.makedirs(Path(config["save_dir"]), exist_ok=True)

    # adjust units.
    config["integrator_config"]["timestep"] *= units.fs
    if config["integrator"] in ['NoseHoover', 'NoseHooverChain']:
        config["integrator_config"]["temperature"] *= units.kB

    # set up simulator.
    integrator = getattr(md_integrator, config["integrator"])(
        atoms, **config["integrator_config"])

    restart = config["read_velocity"]

    simulator = Simulator(atoms, integrator, config["T_init"],
                          restart=restart,
                          start_time=simulated_time,
                          save_dir=Path(config["save_dir"]),
                          save_frequency=config["save_freq"])

    # run simulation.
    start_time = time.time()
    max_std_threshold = config["force_threshold"]
    max_mu_threshold = config["fermi_threshold"]
    early_stop, step = simulator.run(config["steps"] - simulated_step, max_std_threshold, max_mu_threshold)
    elapsed = time.time() - start_time
