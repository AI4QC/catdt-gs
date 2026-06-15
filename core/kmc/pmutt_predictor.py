#!/usr/bin/env python
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'deps', 'pMuTT'))

"""
pMuTT Predictor - A simplified wrapper for multiscale thermochemistry calculations
"""

import json
import tempfile
from typing import Dict, List, Optional, Union, Tuple
from dataclasses import dataclass, field
import numpy as np
import pandas as pd

from pmutt import constants as c
from pmutt.statmech import StatMech, presets
from pmutt.empirical.nasa import Nasa
from pmutt.empirical.references import Reference, References
from pmutt.reaction import Reaction, Reactions
from pmutt.mixture.cov import PiecewiseCovEffect
from pmutt.io import chemkin, cantera as pmutt_cantera
from pmutt.cantera.phase import IdealGas, StoichSolid


# =============================================================================
# Data Classes for Input Specification
# =============================================================================

@dataclass
class SpeciesData:
    """Species thermodynamic data from DFT calculations."""
    name: str                                  # Species name
    atoms: Optional[object] = None             # ASE Atoms object
    potentialenergy: Optional[float] = None    # DFT energy (eV)
    vib_wavenumbers: List[float] = field(default_factory=list)  # cm^-1
    symmetrynumber: int = 1                    # Symmetry number
    spin: float = 0                            # Electronic spin
    phase: str = 'gas'                         # 'gas' or 'surface'
    elements: Dict[str, int] = field(default_factory=dict)  # {'H': 2, 'O': 1}
    geometry: str = 'nonlinear'                # 'nonlinear', 'linear', 'atom'
    # Optional coverage effects
    coverage_effects: List[Dict] = field(default_factory=list)
    # Optional notes
    notes: str = ''


@dataclass
class ReactionData:
    """Reaction data."""
    reactants: List[str]                       # Reactant names
    products: List[str]                        # Product names
    reactants_stoich: List[float] = field(default_factory=list)
    products_stoich: List[float] = field(default_factory=list)
    transition_state: Optional[str] = None     # TS name


@dataclass
class CoverageEffect:
    """Coverage effect specification."""
    affected_species: str                      # Species affected by coverage
    coverage_species: str                      # Species causing coverage
    intervals: List[float]                     # Coverage intervals [0, 0.5, 1.0]
    slopes: List[float]                        # Slopes for each interval (kcal/mol)


# =============================================================================
# Main Predictor Class
# =============================================================================

class pMuTTPredictor:
    """
    Wrapper class for pMuTT thermochemistry calculations.

    Simplifies the process of:
    1. Converting DFT data to NASA polynomials
    2. Applying coverage effects
    3. Calculating reaction thermodynamics
    4. Exporting to Chemkin/Cantera formats

    Example usage:
    -------------
    ```python
    from ase.build import molecule

    predictor = pMuTTPredictor()

    # Add gas species from DFT data
    predictor.add_species(
        name='H2O',
        atoms=molecule('H2O'),
        potentialenergy=-14.2209,  # eV
        vib_wavenumbers=[3825.434, 3710.264, 1582.432],
        phase='gas'
    )

    # Generate NASA polynomials
    nasa_species = predictor.generate_nasa_polynomials()

    # Export to Chemkin
    predictor.export_to_chemkin('thermo.dat')

    # Export to Cantera
    predictor.export_to_cantera('mechanism.yaml')
    ```
    """

    def __init__(self, work_dir: Optional[str] = None, verbose: bool = True):
        """
        Initialize pMuTT predictor.

        Args:
            work_dir: Working directory for output files
            verbose: Print detailed information
        """
        self.work_dir = work_dir or tempfile.mkdtemp(prefix='pmutt_')
        self.verbose = verbose

        # Storage for species and reactions
        self.species_data: Dict[str, SpeciesData] = {}
        self.statmech_models: Dict[str, StatMech] = {}
        self.nasa_models: Dict[str, Nasa] = {}
        self.reactions: List[Reaction] = []

        # Reference species for DFT corrections
        self.references: Optional[References] = None

        # Default settings for NASA polynomial generation
        self.T_low: float = 200.0    # K
        self.T_high: float = 3500.0  # K
        self.T_mid: float = 1000.0   # K

        os.makedirs(self.work_dir, exist_ok=True)

        if self.verbose:
            print(f"pMuTT Predictor initialized:")
            print(f"  Work directory: {self.work_dir}")

    def _log(self, message: str):
        """Print log message if verbose."""
        if self.verbose:
            print(message)

    # =========================================================================
    # Species Management
    # =========================================================================

    def add_species(
        self,
        name: str,
        atoms: Optional[object] = None,
        potentialenergy: Optional[float] = None,
        vib_wavenumbers: Optional[List[float]] = None,
        symmetrynumber: int = 1,
        spin: float = 0,
        phase: str = 'gas',
        elements: Optional[Dict[str, int]] = None,
        geometry: str = 'nonlinear',
        **kwargs
    ):
        """
        Add a species to the database.

        Args:
            name: Species name (e.g., 'H2O', 'CO', 'CO*')
            atoms: ASE Atoms object with structure
            potentialenergy: DFT total energy in eV
            vib_wavenumbers: Vibrational frequencies in cm^-1
            symmetrynumber: Symmetry number (1 for no symmetry)
            spin: Electronic spin
            phase: 'gas' or 'surface'
            elements: Element composition {'H': 2, 'O': 1}
            geometry: 'nonlinear', 'linear', or 'atom'
            **kwargs: Additional parameters
        """
        # Auto-detect elements from atoms if not provided
        if elements is None and atoms is not None:
            try:
                from collections import Counter
                symbols = atoms.get_chemical_symbols()
                elements = dict(Counter(symbols))
            except:
                pass

        # Auto-detect geometry
        if atoms is not None and len(atoms) == 1:
            geometry = 'atom'
        elif atoms is not None and len(atoms) == 2:
            geometry = 'linear'

        species_data = SpeciesData(
            name=name,
            atoms=atoms,
            potentialenergy=potentialenergy,
            vib_wavenumbers=vib_wavenumbers or [],
            symmetrynumber=symmetrynumber,
            spin=spin,
            phase=phase,
            elements=elements or {},
            geometry=geometry,
            **kwargs
        )

        self.species_data[name] = species_data
        self._log(f"Added species: {name} ({phase}, {len(vib_wavenumbers or [])} vibrations)")

    def add_species_from_vasp(
        self,
        name: str,
        outcar_path: str,
        phase: str = 'gas',
        **kwargs
    ):
        """
        Add species from VASP OUTCAR file.

        Args:
            name: Species name
            outcar_path: Path to VASP OUTCAR file
            phase: 'gas' or 'surface'
            **kwargs: Additional parameters
        """
        from pmutt.io.vasp import get_atoms, get_vib_wavenumbers

        atoms = get_atoms(outcar_path)
        vib_wavenumbers = get_vib_wavenumbers(outcar_path)
        potentialenergy = atoms.get_potential_energy()

        self.add_species(
            name=name,
            atoms=atoms,
            potentialenergy=potentialenergy,
            vib_wavenumbers=vib_wavenumbers,
            phase=phase,
            **kwargs
        )

    def add_species_from_gaussian(
        self,
        name: str,
        log_path: str,
        phase: str = 'gas',
        **kwargs
    ):
        """
        Add species from Gaussian log file.

        Args:
            name: Species name
            log_path: Path to Gaussian log file
            phase: 'gas' or 'surface'
            **kwargs: Additional parameters
        """
        from pmutt.io.gaussian import get_atoms, get_vib_wavenumbers

        atoms = get_atoms(log_path)
        vib_wavenumbers = get_vib_wavenumbers(log_path)
        potentialenergy = atoms.get_potential_energy()

        self.add_species(
            name=name,
            atoms=atoms,
            potentialenergy=potentialenergy,
            vib_wavenumbers=vib_wavenumbers,
            phase=phase,
            **kwargs
        )

    def set_reference_species(self, references: List[Dict], T_ref: float = 298.15):
        """
        Set reference species for DFT energy corrections.

        Args:
            references: List of reference species dicts
                Example: [
                    {'name': 'H2', 'elements': {'H': 2}, 'HoRT_ref': 0.0},
                    {'name': 'H2O', 'elements': {'H': 2, 'O': 1}, 'HoRT_ref': -96.6}
                ]
            T_ref: Reference temperature (K), default 298.15

        Note:
            For simple usage, you may not need reference species.
            NASA polynomials can be generated directly from StatMech models.
        """
        ref_objects = []
        for ref_data in references:
            # Add T_ref if not present
            if 'T_ref' not in ref_data:
                ref_data['T_ref'] = T_ref
            ref = Reference(**ref_data)
            ref_objects.append(ref)

        self.references = References(ref_objects)
        self._log(f"Set {len(ref_objects)} reference species for DFT corrections")

    # =========================================================================
    # Coverage Effects
    # =========================================================================

    def add_coverage_effect(
        self,
        affected_species: str,
        coverage_species: str,
        intervals: List[float],
        slopes: List[float]
    ):
        """
        Add coverage effect to a species.

        Args:
            affected_species: Name of species affected by coverage
            coverage_species: Name of species causing coverage
            intervals: Coverage intervals [0, 0.3, 0.7, 1.0]
            slopes: Slopes for each interval in kcal/mol [-8.0, -5.0, -2.0]

        Example:
            # O coverage affects CO adsorption energy
            predictor.add_coverage_effect(
                affected_species='CO*',
                coverage_species='O*',
                intervals=[0.0, 0.3, 0.7, 1.0],
                slopes=[-8.0, -5.0, -2.0]
            )
        """
        if affected_species not in self.species_data:
            raise ValueError(f"Species '{affected_species}' not found. Add it first using add_species().")

        cov_effect = PiecewiseCovEffect(
            name_i=affected_species,
            name_j=coverage_species,
            intervals=intervals,
            slopes=slopes
        )

        # Add to species data
        self.species_data[affected_species].coverage_effects.append({
            'coverage_species': coverage_species,
            'intervals': intervals,
            'slopes': slopes,
            'model': cov_effect
        })

        self._log(f"Added coverage effect: {coverage_species} affects {affected_species}")

    # =========================================================================
    # StatMech Model Generation
    # =========================================================================

    def generate_statmech_models(self):
        """
        Generate statistical mechanics models for all species.

        Returns:
            Dict of StatMech models {name: model}
        """
        self._log("\nGenerating statistical mechanics models...")

        for name, data in self.species_data.items():
            # Choose preset based on phase
            if data.phase == 'gas':
                if data.geometry == 'atom':
                    preset = presets['electronic']
                else:
                    preset = presets['idealgas']
            else:  # surface
                preset = presets['harmonic']

            # Create StatMech model
            statmech_kwargs = {
                'name': name,
                'potentialenergy': data.potentialenergy,
                'symmetrynumber': data.symmetrynumber,
                'spin': data.spin,
                **preset
            }

            # Add atoms if available
            if data.atoms is not None:
                statmech_kwargs['atoms'] = data.atoms

            # Add vibrations
            if data.vib_wavenumbers:
                statmech_kwargs['vib_wavenumbers'] = data.vib_wavenumbers

            # Create model
            model = StatMech(**statmech_kwargs)

            # Apply coverage effects if any
            if data.coverage_effects:
                misc_models = [eff['model'] for eff in data.coverage_effects]
                model.misc_models = misc_models

            self.statmech_models[name] = model
            self._log(f"  Created StatMech model for {name}")

        return self.statmech_models

    # =========================================================================
    # NASA Polynomial Generation
    # =========================================================================

    def generate_nasa_polynomials(
        self,
        T_low: Optional[float] = None,
        T_high: Optional[float] = None,
        T_mid: Optional[float] = None
    ) -> Dict[str, Nasa]:
        """
        Generate NASA polynomials from statistical mechanics models.

        Args:
            T_low: Lower temperature limit (K), default 200
            T_high: Upper temperature limit (K), default 3500
            T_mid: Mid temperature (K), default 1000

        Returns:
            Dict of NASA models {name: model}
        """
        if not self.statmech_models:
            self.generate_statmech_models()

        T_low = T_low or self.T_low
        T_high = T_high or self.T_high
        T_mid = T_mid or self.T_mid

        self._log(f"\nGenerating NASA polynomials (T: {T_low}-{T_mid}-{T_high} K)...")

        for name, statmech_model in self.statmech_models.items():
            data = self.species_data[name]

            # Generate NASA polynomial using from_model (replaces deprecated from_statmech)
            nasa = Nasa.from_model(
                name=name,
                model=statmech_model,
                T_low=T_low,
                T_high=T_high,
                T_mid=T_mid,
                elements=data.elements,
                phase=data.phase,
                references=self.references
            )

            self.nasa_models[name] = nasa
            self._log(f"  Generated NASA polynomial for {name}")

        return self.nasa_models

    # =========================================================================
    # Reaction Management
    # =========================================================================

    def add_reaction(
        self,
        reactants: List[str],
        products: List[str],
        reactants_stoich: Optional[List[float]] = None,
        products_stoich: Optional[List[float]] = None,
        transition_state: Optional[str] = None
    ):
        """
        Add a reaction.

        Args:
            reactants: List of reactant species names
            products: List of product species names
            reactants_stoich: Stoichiometric coefficients (default: all 1)
            products_stoich: Stoichiometric coefficients (default: all 1)
            transition_state: Name of transition state species

        Example:
            predictor.add_reaction(
                reactants=['CO*', 'O*'],
                products=['CO2', '*', '*'],
                transition_state='CO-O*'
            )
        """
        # Default stoichiometry
        if reactants_stoich is None:
            reactants_stoich = [1.0] * len(reactants)
        if products_stoich is None:
            products_stoich = [1.0] * len(products)

        # Get species objects
        reactant_species = []
        for r_name in reactants:
            if r_name in self.nasa_models:
                reactant_species.append(self.nasa_models[r_name])
            elif r_name in self.statmech_models:
                reactant_species.append(self.statmech_models[r_name])
            else:
                raise ValueError(f"Reactant '{r_name}' not found in models")

        product_species = []
        for p_name in products:
            if p_name in self.nasa_models:
                product_species.append(self.nasa_models[p_name])
            elif p_name in self.statmech_models:
                product_species.append(self.statmech_models[p_name])
            else:
                raise ValueError(f"Product '{p_name}' not found in models")

        # Transition state
        ts_species = None
        if transition_state:
            if transition_state in self.nasa_models:
                ts_species = self.nasa_models[transition_state]
            elif transition_state in self.statmech_models:
                ts_species = self.statmech_models[transition_state]

        # Create reaction
        reaction = Reaction(
            reactants=reactant_species,
            reactants_stoich=reactants_stoich,
            products=product_species,
            products_stoich=products_stoich,
            transition_state=ts_species
        )

        self.reactions.append(reaction)
        self._log(f"Added reaction: {' + '.join(reactants)} -> {' + '.join(products)}")

    def calculate_reaction_properties(
        self,
        temperature: float = 298.15,
        units: str = 'kJ/mol'
    ) -> pd.DataFrame:
        """
        Calculate thermodynamic properties for all reactions.

        Args:
            temperature: Temperature in K
            units: Energy units ('kJ/mol', 'eV', 'kcal/mol')

        Returns:
            DataFrame with reaction properties
        """
        results = []

        for i, rxn in enumerate(self.reactions):
            result = {
                'reaction_index': i,
                'temperature': temperature,
                'delta_H': rxn.get_delta_H(units=units, T=temperature),
                'delta_S': rxn.get_delta_S(units=f'{units}/K', T=temperature),
                'delta_G': rxn.get_delta_G(units=units, T=temperature),
            }

            # Activation energies if TS is available
            if rxn.transition_state is not None:
                result['Ea_forward'] = rxn.get_E_act(rev=False, units=units, T=temperature)
                result['Ea_reverse'] = rxn.get_E_act(rev=True, units=units, T=temperature)

            results.append(result)

        return pd.DataFrame(results)

    # =========================================================================
    # Export Functions
    # =========================================================================

    def export_to_chemkin(
        self,
        filename: str = 'thermo.dat',
        species_list: Optional[List[str]] = None
    ):
        """
        Export NASA polynomials to Chemkin THERMO.DAT format.

        Args:
            filename: Output filename
            species_list: List of species to export (None = all)
        """
        if not self.nasa_models:
            self.generate_nasa_polynomials()

        # Select species
        if species_list is None:
            nasa_list = list(self.nasa_models.values())
        else:
            nasa_list = [self.nasa_models[name] for name in species_list]

        # Export
        output_path = os.path.join(self.work_dir, filename)
        chemkin.write_gas(nasa_list, filename=output_path)

        self._log(f"\nExported {len(nasa_list)} species to Chemkin format: {output_path}")
        return output_path

    def export_to_cantera(
        self,
        filename: str = 'mechanism.yaml',
        gas_species: Optional[List[str]] = None,
        surface_species: Optional[List[str]] = None
    ):
        """
        Export to Cantera YAML format.

        Args:
            filename: Output filename
            gas_species: List of gas species names (None = auto-detect)
            surface_species: List of surface species names (None = auto-detect)
        """
        if not self.nasa_models:
            self.generate_nasa_polynomials()

        # Auto-detect phases if not specified
        if gas_species is None:
            gas_species = [name for name, data in self.species_data.items()
                          if data.phase == 'gas']
        if surface_species is None:
            surface_species = [name for name, data in self.species_data.items()
                              if data.phase == 'surface']

        phases = []

        # Gas phase
        if gas_species:
            gas_nasa = [self.nasa_models[name] for name in gas_species]
            gas_phase = IdealGas(
                name='gas',
                species=gas_nasa,
                transport='Mix',
                kinetics='GasKinetics'
            )
            phases.append(gas_phase)

        # Surface phase
        if surface_species:
            surf_nasa = [self.nasa_models[name] for name in surface_species]
            surf_phase = StoichSolid(
                name='surface',
                species=surf_nasa,
                density=2.7e-9  # mol/cm^2, typical for metal surfaces
            )
            phases.append(surf_phase)

        # Convert to CTI string
        cti_strings = []
        for phase in phases:
            cti_str = pmutt_cantera.obj_to_cti(phase)
            cti_strings.append(cti_str)

        # Write to file
        output_path = os.path.join(self.work_dir, filename)
        with open(output_path, 'w') as f:
            f.write('\n\n'.join(cti_strings))

        self._log(f"\nExported to Cantera format: {output_path}")
        self._log(f"  Gas species: {len(gas_species)}")
        self._log(f"  Surface species: {len(surface_species)}")

        return output_path

    def export_to_excel(
        self,
        filename: str = 'thermodynamics.xlsx',
        temperature_range: Optional[List[float]] = None
    ):
        """
        Export thermodynamic data to Excel file.

        Args:
            filename: Output filename
            temperature_range: List of temperatures to evaluate (K)
        """
        if not self.nasa_models:
            self.generate_nasa_polynomials()

        if temperature_range is None:
            temperature_range = [298.15, 500, 1000, 1500, 2000]

        # Prepare data
        data = []
        for name, nasa in self.nasa_models.items():
            for T in temperature_range:
                row = {
                    'Species': name,
                    'Phase': self.species_data[name].phase,
                    'Temperature (K)': T,
                    'Cp (J/mol/K)': nasa.get_Cp(units='J/mol/K', T=T),
                    'H (kJ/mol)': nasa.get_H(units='kJ/mol', T=T),
                    'S (J/mol/K)': nasa.get_S(units='J/mol/K', T=T),
                    'G (kJ/mol)': nasa.get_G(units='kJ/mol', T=T),
                }
                data.append(row)

        df = pd.DataFrame(data)

        # Write to Excel
        output_path = os.path.join(self.work_dir, filename)
        df.to_excel(output_path, index=False)

        self._log(f"\nExported thermodynamic data to Excel: {output_path}")
        return output_path

    # =========================================================================
    # Calculation Methods
    # =========================================================================

    def calculate_thermo_at_T(
        self,
        species_name: str,
        temperature: float,
        units: str = 'kJ/mol'
    ) -> Dict[str, float]:
        """
        Calculate thermodynamic properties for a species at given temperature.

        Args:
            species_name: Name of species
            temperature: Temperature in K
            units: Energy units

        Returns:
            Dict with Cp, H, S, G
        """
        if species_name in self.nasa_models:
            model = self.nasa_models[species_name]
        elif species_name in self.statmech_models:
            model = self.statmech_models[species_name]
        else:
            raise ValueError(f"Species '{species_name}' not found")

        return {
            'Cp': model.get_Cp(T=temperature, units=f'{units}/K'),
            'H': model.get_H(T=temperature, units=units),
            'S': model.get_S(T=temperature, units=f'{units}/K'),
            'G': model.get_G(T=temperature, units=units),
        }

    def get_species_summary(self) -> pd.DataFrame:
        """Get summary of all species."""
        data = []
        for name, spec_data in self.species_data.items():
            row = {
                'Name': name,
                'Phase': spec_data.phase,
                'Elements': str(spec_data.elements),
                'Energy (eV)': spec_data.potentialenergy,
                'N_vibrations': len(spec_data.vib_wavenumbers),
                'Symmetry': spec_data.symmetrynumber,
                'Coverage_effects': len(spec_data.coverage_effects),
                'Has_NASA': name in self.nasa_models,
            }
            data.append(row)

        return pd.DataFrame(data)

    def summary(self) -> str:
        """Generate a summary of the predictor state."""
        lines = ["=" * 60]
        lines.append("pMuTT Predictor Summary")
        lines.append("=" * 60)
        lines.append(f"Work directory: {self.work_dir}")
        lines.append(f"")
        lines.append(f"Species: {len(self.species_data)}")

        # Count by phase
        gas_count = sum(1 for d in self.species_data.values() if d.phase == 'gas')
        surf_count = sum(1 for d in self.species_data.values() if d.phase == 'surface')
        lines.append(f"  Gas: {gas_count}")
        lines.append(f"  Surface: {surf_count}")

        lines.append(f"")
        lines.append(f"Models:")
        lines.append(f"  StatMech: {len(self.statmech_models)}")
        lines.append(f"  NASA: {len(self.nasa_models)}")
        lines.append(f"")
        lines.append(f"Reactions: {len(self.reactions)}")

        if self.references:
            lines.append(f"Reference species: {len(self.references.references)}")

        lines.append("=" * 60)
        return "\n".join(lines)


# =============================================================================
# Convenience Functions
# =============================================================================

def quick_nasa_from_dft(
    name: str,
    atoms: object,
    potentialenergy: float,
    vib_wavenumbers: List[float],
    phase: str = 'gas',
    **kwargs
) -> Nasa:
    """
    Quickly generate NASA polynomial from DFT data.

    Args:
        name: Species name
        atoms: ASE Atoms object
        potentialenergy: DFT energy (eV)
        vib_wavenumbers: Vibrational frequencies (cm^-1)
        phase: 'gas' or 'surface'
        **kwargs: Additional parameters

    Returns:
        Nasa object
    """
    predictor = pMuTTPredictor(verbose=False)
    predictor.add_species(
        name=name,
        atoms=atoms,
        potentialenergy=potentialenergy,
        vib_wavenumbers=vib_wavenumbers,
        phase=phase,
        **kwargs
    )
    predictor.generate_statmech_models()
    nasa_models = predictor.generate_nasa_polynomials()
    return nasa_models[name]


# =============================================================================
# Main entry point for testing
# =============================================================================

if __name__ == '__main__':
    print("pMuTT Predictor - Testing Example")
    print("=" * 60)

    from ase.build import molecule

    # Create predictor
    predictor = pMuTTPredictor(work_dir='output/pmutt_test')

    # =========================================================================
    # Test 1: Add gas species
    # =========================================================================
    print("\nTest 1: Adding gas species")
    print("-" * 60)

    # H2
    predictor.add_species(
        name='H2',
        atoms=molecule('H2'),
        potentialenergy=-6.77,
        vib_wavenumbers=[4342],
        symmetrynumber=2,
        phase='gas'
    )

    # O2
    predictor.add_species(
        name='O2',
        atoms=molecule('O2'),
        potentialenergy=-9.86,
        vib_wavenumbers=[1580],
        symmetrynumber=2,
        spin=1,
        phase='gas'
    )

    # H2O
    predictor.add_species(
        name='H2O',
        atoms=molecule('H2O'),
        potentialenergy=-14.22,
        vib_wavenumbers=[3825.434, 3710.264, 1582.432],
        symmetrynumber=2,
        phase='gas'
    )

    # =========================================================================
    # Test 2: Generate NASA polynomials (no references needed for basic usage)
    # =========================================================================
    print("\nTest 2: Generating NASA polynomials")
    print("-" * 60)

    nasa_models = predictor.generate_nasa_polynomials()

    # Calculate properties at 500 K
    print("\nThermodynamic properties at 500 K:")
    for name in ['H2', 'O2', 'H2O']:
        props = predictor.calculate_thermo_at_T(name, 500.0, units='kJ/mol')
        print(f"\n{name}:")
        print(f"  Cp = {props['Cp']:.2f} kJ/mol/K")
        print(f"  H  = {props['H']:.2f} kJ/mol")
        print(f"  S  = {props['S']:.2f} kJ/mol/K")
        print(f"  G  = {props['G']:.2f} kJ/mol")

    # =========================================================================
    # Test 3: Add reaction and calculate properties
    # =========================================================================
    print("\nTest 3: Reaction thermodynamics")
    print("-" * 60)

    predictor.add_reaction(
        reactants=['H2', 'O2'],
        products=['H2O'],
        reactants_stoich=[1.0, 0.5],
        products_stoich=[1.0]
    )

    rxn_props = predictor.calculate_reaction_properties(temperature=500.0)
    print("\nH2 + 0.5 O2 -> H2O at 500 K:")
    print(rxn_props.to_string(index=False))

    # =========================================================================
    # Test 4: Export to different formats
    # =========================================================================
    print("\nTest 4: Exporting to different formats")
    print("-" * 60)

    # Chemkin
    chemkin_file = predictor.export_to_chemkin()
    print(f"\nChemkin file: {chemkin_file}")

    # Cantera
    cantera_file = predictor.export_to_cantera()
    print(f"Cantera file: {cantera_file}")

    # Excel
    excel_file = predictor.export_to_excel()
    print(f"Excel file: {excel_file}")

    # =========================================================================
    # Test 5: Summary
    # =========================================================================
    print("\n" + predictor.summary())

    print("\n" + "=" * 60)
    print("Test completed successfully!")
    print("=" * 60)
