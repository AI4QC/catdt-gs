#!/usr/bin/env python
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'deps', 'catmap'))

"""
CatMAP Predictor - A wrapper for microkinetic modeling using CatMAP
"""

import tempfile
import shutil
from typing import Any, Dict, List, Optional, Union, Tuple
from dataclasses import dataclass, field
import numpy as np

# =============================================================================
# Data Classes for Input Specification
# =============================================================================

@dataclass
class GasSpecies:
    """Gas phase species definition."""
    name: str                          # e.g., 'CO', 'O2', 'H2O'
    pressure: float = 1.0              # Partial pressure (bar)
    formation_energy: float = 0.0      # Formation energy (eV), relative to references
    frequencies: List[float] = field(default_factory=list)  # Vibrational frequencies (cm^-1)
    # Optional ideal-gas thermodynamic parameters for CatMAP's ideal_gas mode.
    # None for any entry → CatMAP falls back to its built-in lookup table;
    # unknown molecules need these supplied to compute TΔS(T).
    symmetry_number: Optional[int] = None
    geometry: Optional[str] = None     # 'linear' or 'nonlinear' or 'monatomic'
    spin: Optional[float] = None       # e.g. 0 for H2, 1.0 for O2 (triplet)
    # Relaxed gas-phase geometry (used by CatMAP to compute moments of inertia
    # for the rotational partition function). If unset, CatMAP falls back to
    # `ase.build.molecule(name)` which fails for molecules outside ASE's g2
    # collection (e.g. propane).
    atoms_geometry: Optional[Any] = None  # ASE Atoms

@dataclass
class AdsorbateSpecies:
    """Adsorbed species definition."""
    name: str                          # e.g., 'CO', 'O', 'OH'
    site: str = '111'                  # Adsorption site
    formation_energies: Dict[str, float] = field(default_factory=dict)  # {surface: energy}
    frequencies: List[float] = field(default_factory=list)

@dataclass
class TransitionState:
    """Transition state definition."""
    name: str                          # e.g., 'O-CO', 'H-OH'
    site: str = '111'
    formation_energies: Dict[str, float] = field(default_factory=dict)  # {surface: energy}
    frequencies: List[float] = field(default_factory=list)
    scaling_mode: str = 'initial_state'  # 'initial_state', 'final_state', 'BEP', or list

@dataclass
class ReactionStep:
    """Elementary reaction step."""
    expression: str                    # e.g., '*_s + CO_g -> CO*'
    # Supports CatMAP syntax:
    # - Adsorption: '*_s + CO_g -> CO*'
    # - Dissociation: '2*_s + O2_g <-> O-O* + *_s -> 2O*'
    # - Surface reaction: 'CO* + O* <-> O-CO* + * -> CO2_g + 2*'
    # - Electrochemical: 'H* + pe_g -> H2_g + *_s' (pe = proton-electron pair)


# =============================================================================
# Main Predictor Class
# =============================================================================

class CatMAPPredictor:
    """
    Wrapper class for CatMAP microkinetic modeling.

    Simplifies the process of setting up and running microkinetic models
    for arbitrary reaction systems.

    Example usage:
    -------------
    ```python
    predictor = CatMAPPredictor()

    # Define reaction mechanism
    predictor.add_reaction('*_s + CO_g -> CO*')
    predictor.add_reaction('2*_s + O2_g <-> O-O* + *_s -> 2O*')
    predictor.add_reaction('CO* + O* <-> O-CO* + * -> CO2_g + 2*')

    # Add energy data
    predictor.add_gas('CO', formation_energy=0.0, frequencies=[2170])
    predictor.add_gas('O2', formation_energy=0.0, frequencies=[1580])
    predictor.add_gas('CO2', formation_energy=-2.45, frequencies=[1333, 2349, 667, 667])

    predictor.add_adsorbate('CO', {'Pt': 1.7, 'Pd': 1.55, 'Cu': 2.58})
    predictor.add_adsorbate('O', {'Pt': 1.62, 'Pd': 1.55, 'Cu': 1.07})
    predictor.add_transition_state('O-CO', {'Pt': 4.04, 'Pd': 4.2, 'Cu': 4.18})

    # Set conditions and run
    predictor.set_conditions(temperature=500, pressures={'CO_g': 1.0, 'O2_g': 0.33})
    result = predictor.run()
    ```
    """

    # Available thermodynamic modes
    GAS_THERMO_MODES = [
        'shomate_gas',      # Shomate polynomials (most accurate)
        'ideal_gas',        # Ideal gas with statistical mechanics
        'zero_point_gas',   # Only zero-point corrections
        'fixed_entropy_gas', # Fixed entropy approximation
        'frozen_gas',       # No thermal contributions
    ]

    ADSORBATE_THERMO_MODES = [
        'harmonic_adsorbate',   # Harmonic oscillator
        'hindered_adsorbate',   # Hindered translator/rotor
        'zero_point_adsorbate', # Only zero-point corrections
        'frozen_adsorbate',     # No thermal contributions
    ]

    ELECTROCHEMICAL_THERMO_MODES = [
        'simple_electrochemical',    # Basic CHE model
        'hbond_electrochemical',     # With H-bonding corrections
        'hbond_with_estimates_electrochemical',
        'field_electrochemical',     # With field effects
        'local_field_electrochemical',
    ]

    # Available output variables
    OUTPUT_VARIABLES = [
        # Solver level outputs
        'coverage',           # Surface coverages
        'rate',               # Reaction rates
        'production_rate',    # Production rates of species
        'consumption_rate',   # Consumption rates
        'turnover_frequency', # TOF
        'selectivity',        # Product selectivity
        'carbon_selectivity', # Carbon-based selectivity
        'rate_control',       # Degree of rate control
        'selectivity_control',# Degree of selectivity control
        'rxn_direction',      # Net reaction direction
        'directional_rates',  # Forward/reverse rates
        'rxn_order',          # Apparent reaction orders
        'apparent_activation_energy',
        # Scaler level outputs
        'rxn_parameter',      # Reaction parameters
        'frequency',          # Frequencies
        'electronic_energy',  # Electronic energies
        'free_energy',        # Free energies
        'zero_point_energy',  # ZPE
        'enthalpy',           # Enthalpies
        'entropy',            # Entropies
    ]

    def __init__(self, catmap_root: Optional[str] = None, work_dir: Optional[str] = None):
        """
        Initialize CatMAP predictor.

        Args:
            catmap_root: Path to CatMAP installation (auto-detected if None)
            work_dir: Working directory for output files
        """
        self.catmap_root = catmap_root
        self.work_dir = work_dir or tempfile.mkdtemp(prefix='catmap_')

        # Storage for model components
        self.reactions: List[str] = []
        self.gases: Dict[str, GasSpecies] = {}
        self.adsorbates: Dict[str, AdsorbateSpecies] = {}
        self.transition_states: Dict[str, TransitionState] = {}
        self.surfaces: List[str] = []

        # Model settings
        self.temperature: float = 500.0  # K
        self.voltage: Optional[float] = None  # V vs SHE (for electrochemistry)
        self.beta: float = 0.5  # Symmetry factor
        self.pH: float = 0.0

        # Descriptor settings
        self.descriptor_names: List[str] = []
        self.descriptor_ranges: List[List[float]] = []
        self.resolution: int = 15

        # Thermodynamic settings
        self.gas_thermo_mode: str = 'ideal_gas'
        self.adsorbate_thermo_mode: str = 'harmonic_adsorbate'
        self.electrochemical_thermo_mode: str = 'simple_electrochemical'

        # Scaling settings
        self.scaling_constraints: Dict[str, Union[str, List]] = {}

        # Interaction settings
        self.adsorbate_interaction_model: Optional[str] = None  # None, 'ideal', 'first_order'
        self.interaction_parameters: Dict[str, Dict] = {}

        # Output settings
        self.output_variables: List[str] = ['coverage', 'rate', 'production_rate']

        # Site definition
        self.site_name: str = 's'
        self.site_type: str = '111'
        self.total_sites: float = 1.0

        # Atomic reservoir settings (for defining references for each element)
        self.atomic_reservoir_list: List[str] = []

        # Solver settings
        self.decimal_precision: int = int(os.getenv("CATDT_CATMAP_DECIMAL_PRECISION", "100"))
        self.tolerance: float = float(os.getenv("CATDT_CATMAP_TOLERANCE", "1e-50"))
        self.max_iterations: int = int(os.getenv("CATDT_CATMAP_MAX_ITERATIONS", "100"))
        use_numbers = os.getenv("CATDT_CATMAP_USE_NUMBERS_SOLVER")
        self.use_numbers_solver: Optional[bool] = None
        if use_numbers is not None:
            self.use_numbers_solver = use_numbers.strip().lower() in {"1", "true", "yes", "on"}
        self.max_bisections: int = int(os.getenv("CATDT_CATMAP_MAX_BISECTIONS", "5"))
        self.max_damping_iterations: Optional[int] = None
        if os.getenv("CATDT_CATMAP_MAX_DAMPING_ITERATIONS") is not None:
            self.max_damping_iterations = int(os.getenv("CATDT_CATMAP_MAX_DAMPING_ITERATIONS", "10"))
        self.max_initial_guesses: Optional[int] = None
        if os.getenv("CATDT_CATMAP_MAX_INITIAL_GUESSES") is not None:
            self.max_initial_guesses = int(os.getenv("CATDT_CATMAP_MAX_INITIAL_GUESSES", "3"))

        # Results storage
        self.model = None
        self.results = None

    # =========================================================================
    # Methods for adding model components
    # =========================================================================

    def add_reaction(self, expression: str):
        """
        Add an elementary reaction step.

        Supports CatMAP syntax:
        - Adsorption: '*_s + CO_g -> CO*' or '*_s + CO_g -> CO_s'
        - Desorption: 'CO* -> *_s + CO_g'
        - Dissociation with TS: '2*_s + O2_g <-> O-O* + *_s -> 2O*'
        - Surface reaction: 'CO* + O* <-> O-CO* + * -> CO2_g + 2*'
        - Electrochemical: 'H* + pe_g -> H2_g + *_s' (pe = proton-electron)

        Args:
            expression: Reaction expression string
        """
        self.reactions.append(expression)

    def add_gas(self, name: str,
                formation_energy: float = 0.0,
                pressure: float = 1.0,
                frequencies: Optional[List[float]] = None,
                symmetry_number: Optional[int] = None,
                geometry: Optional[str] = None,
                spin: Optional[float] = None,
                atoms_geometry: Optional[Any] = None):
        """
        Add a gas phase species.

        Args:
            name: Species name (e.g., 'CO', 'O2', 'H2O'). Trailing '_g' is
                stripped so callers may pass either the bare name or the
                suffixed form; _generate_single_point_mkm always appends '_g'
                once for the CatMAP species key.
            formation_energy: Formation energy in eV (relative to references)
            pressure: Partial pressure in bar
            frequencies: Vibrational frequencies in cm^-1 (required for ideal_gas)
            symmetry_number, geometry, spin: Ideal-gas parameters used by
                CatMAP's ideal_gas thermo mode. See
                catmap.data.parameter_data.ideal_gas_params for the built-in
                table; unknown molecules must supply these explicitly.
        """
        if name.endswith("_g"):
            name = name[:-2]
        self.gases[name] = GasSpecies(
            name=name,
            formation_energy=formation_energy,
            pressure=pressure,
            frequencies=frequencies or [],
            symmetry_number=symmetry_number,
            geometry=geometry,
            spin=spin,
            atoms_geometry=atoms_geometry,
        )

    def add_adsorbate(self, name: str,
                      formation_energies: Dict[str, float],
                      site: str = '111',
                      frequencies: Optional[List[float]] = None):
        """
        Add an adsorbed species.

        Args:
            name: Adsorbate name (e.g., 'CO', 'O', 'OH'). Trailing '_s' is
                stripped; the CatMAP species key is reconstructed as
                '<name>_s' in the generated mkm file.
            formation_energies: Dict mapping surface names to formation energies (eV)
            site: Adsorption site type (default: '111')
            frequencies: Vibrational frequencies in cm^-1
        """
        if name.endswith("_s"):
            name = name[:-2]
        self.adsorbates[name] = AdsorbateSpecies(
            name=name,
            site=site,
            formation_energies=formation_energies,
            frequencies=frequencies or []
        )
        # Auto-add surfaces
        for surf in formation_energies.keys():
            if surf not in self.surfaces:
                self.surfaces.append(surf)

    def add_transition_state(self, name: str,
                             formation_energies: Dict[str, float],
                             site: str = '111',
                             frequencies: Optional[List[float]] = None,
                             scaling_mode: str = 'initial_state'):
        """
        Add a transition state.

        Args:
            name: TS name (e.g., 'O-CO', 'H-OH') - use '-' to indicate bond forming/breaking.
                Trailing '_s' is stripped; the CatMAP species key is reconstructed
                as '<name>_s' in the generated mkm file.
            formation_energies: Dict mapping surface names to formation energies (eV)
            site: Site type
            frequencies: Vibrational frequencies
            scaling_mode: 'initial_state', 'final_state', 'BEP', or explicit scaling list
        """
        if name.endswith("_s"):
            name = name[:-2]
        self.transition_states[name] = TransitionState(
            name=name,
            site=site,
            formation_energies=formation_energies,
            frequencies=frequencies or [],
            scaling_mode=scaling_mode
        )
        # Auto-add surfaces
        for surf in formation_energies.keys():
            if surf not in self.surfaces:
                self.surfaces.append(surf)

    # =========================================================================
    # Methods for setting conditions
    # =========================================================================

    def set_conditions(self,
                       temperature: float = 500.0,
                       pressures: Optional[Dict[str, float]] = None,
                       voltage: Optional[float] = None,
                       pH: float = 0.0,
                       beta: float = 0.5):
        """
        Set reaction conditions.

        Args:
            temperature: Temperature in K
            pressures: Dict of gas partial pressures {species_g: pressure}
            voltage: Electrode potential in V vs SHE (for electrochemistry)
            pH: Solution pH (for electrochemistry)
            beta: Symmetry factor for electrochemical reactions
        """
        self.temperature = temperature
        self.voltage = voltage
        self.pH = pH
        self.beta = beta

        if pressures:
            for species, p in pressures.items():
                # Remove '_g' suffix if present
                name = species.replace('_g', '')
                if name in self.gases:
                    self.gases[name].pressure = p

    def set_atomic_reservoirs(self, reservoir_list: List[str]):
        """
        Set the atomic reservoir list for CatMAP.

        This tells CatMAP which gas-phase species to use as references for
        each element in the reaction mechanism.

        Args:
            reservoir_list: List of gas species (e.g., ['H2_g', 'CO_g', 'H2O_g'])
        """
        self.atomic_reservoir_list = reservoir_list

    def set_descriptors(self,
                        names: List[str],
                        ranges: List[List[float]],
                        resolution: int = 15):
        """
        Set descriptor space for volcano plots.

        Args:
            names: Descriptor names (must match adsorbate names, e.g., ['O_s', 'CO_s'])
            ranges: Descriptor ranges [[min1, max1], [min2, max2]]
            resolution: Grid resolution
        """
        self.descriptor_names = names
        self.descriptor_ranges = ranges
        self.resolution = resolution

    def set_thermo_modes(self,
                         gas_mode: str = 'ideal_gas',
                         adsorbate_mode: str = 'harmonic_adsorbate',
                         electrochemical_mode: str = 'simple_electrochemical'):
        """
        Set thermodynamic correction modes.

        Args:
            gas_mode: One of GAS_THERMO_MODES
            adsorbate_mode: One of ADSORBATE_THERMO_MODES
            electrochemical_mode: One of ELECTROCHEMICAL_THERMO_MODES
        """
        if gas_mode not in self.GAS_THERMO_MODES:
            raise ValueError(f"Invalid gas_mode. Choose from: {self.GAS_THERMO_MODES}")
        if adsorbate_mode not in self.ADSORBATE_THERMO_MODES:
            raise ValueError(f"Invalid adsorbate_mode. Choose from: {self.ADSORBATE_THERMO_MODES}")

        self.gas_thermo_mode = gas_mode
        self.adsorbate_thermo_mode = adsorbate_mode
        self.electrochemical_thermo_mode = electrochemical_mode

    def set_scaling_constraints(self, constraints: Dict[str, Union[str, List]]):
        """
        Set scaling constraints for adsorbates and transition states.

        Args:
            constraints: Dict mapping species to constraints
                - For adsorbates: [slope1, slope2, intercept] or ['+', 0, None]
                  '+' means fit slope, 0 means fixed at 0, None means fit intercept
                - For TS: 'initial_state', 'final_state', 'BEP', or explicit list

        Example:
            predictor.set_scaling_constraints({
                'O_s': ['+', 0, None],      # Scales only with descriptor 1
                'CO_s': [0, '+', None],     # Scales only with descriptor 2
                'O-CO_s': 'initial_state',  # TS scales with IS
            })
        """
        self.scaling_constraints = constraints

    def set_interactions(self,
                         model: str = 'first_order',
                         self_interactions: Optional[Dict[str, List[float]]] = None,
                         cross_interactions: Optional[Dict[str, Dict[str, List[float]]]] = None):
        """
        Set adsorbate-adsorbate interaction parameters.

        Args:
            model: 'first_order' or None
            self_interactions: {species: [param_per_surface]}
            cross_interactions: {species1: {species2: [param_per_surface]}}
        """
        self.adsorbate_interaction_model = model
        if self_interactions:
            for sp, params in self_interactions.items():
                if sp not in self.interaction_parameters:
                    self.interaction_parameters[sp] = {}
                self.interaction_parameters[sp]['self_interaction_parameter'] = params
        if cross_interactions:
            for sp1, cross_dict in cross_interactions.items():
                if sp1 not in self.interaction_parameters:
                    self.interaction_parameters[sp1] = {}
                self.interaction_parameters[sp1]['cross_interaction_parameters'] = cross_dict

    # =========================================================================
    # Input file generation
    # =========================================================================

    def _generate_energies_file(self, filepath: str):
        """Generate the energies.txt input file.

        CatMAP expects `frequencies` in **eV** (passed directly to
        ase.thermochemistry.IdealGasThermo / HarmonicThermo, which consume
        vib_energies in eV). Callers should convert from cm^-1 (0.0001240 eV
        per cm^-1) before calling add_gas/add_adsorbate/add_transition_state.
        """
        lines = ['surface_name\tsite_name\tspecies_name\tformation_energy\tbulk_structure\tfrequencies\tother_parameters\treference']

        # Add gas species
        for name, gas in self.gases.items():
            freq_str = str(list(gas.frequencies)) if gas.frequencies else '[]'
            lines.append(f"None\tgas\t{name}\t{gas.formation_energy}\tNone\t{freq_str}\t[]\tCatMAPPredictor")

        # Add adsorbates
        for name, ads in self.adsorbates.items():
            freq_str = str(ads.frequencies) if ads.frequencies else '[]'
            for surf, energy in ads.formation_energies.items():
                lines.append(f"{surf}\t{ads.site}\t{name}\t{energy}\tfcc\t{freq_str}\t[]\tCatMAPPredictor")

        # Add transition states
        for name, ts in self.transition_states.items():
            freq_str = str(ts.frequencies) if ts.frequencies else '[]'
            for surf, energy in ts.formation_energies.items():
                lines.append(f"{surf}\t{ts.site}\t{name}\t{energy}\tfcc\t{freq_str}\t[]\tCatMAPPredictor")

        with open(filepath, 'w') as f:
            f.write('\n'.join(lines))

    def _generate_mkm_file(self, filepath: str, energies_file: str):
        """Generate the .mkm setup file."""
        lines = []

        # Reaction expressions
        lines.append("rxn_expressions = [")
        for rxn in self.reactions:
            lines.append(f"    '{rxn}',")
        lines.append("]")
        lines.append("")

        # Surface names
        lines.append(f"surface_names = {self.surfaces}")
        lines.append("")

        # Descriptor names and ranges
        if self.descriptor_names:
            lines.append(f"descriptor_names = {self.descriptor_names}")
            lines.append(f"descriptor_ranges = {self.descriptor_ranges}")
            lines.append(f"resolution = {self.resolution}")
        else:
            # Auto-detect descriptors from adsorbates (use first two)
            ads_names = list(self.adsorbates.keys())[:2]
            desc_names = [f"{name}_s" for name in ads_names]
            lines.append(f"descriptor_names = {desc_names}")
            # Auto-determine ranges from formation energies
            ranges = []
            for name in ads_names:
                energies = list(self.adsorbates[name].formation_energies.values())
                ranges.append([min(energies) - 1, max(energies) + 1])
            lines.append(f"descriptor_ranges = {ranges}")
            lines.append(f"resolution = {self.resolution}")
        lines.append("")

        # Temperature
        lines.append(f"temperature = {self.temperature}")
        lines.append("")

        # Electrochemical settings
        if self.voltage is not None:
            lines.append(f"voltage = {self.voltage}")
            lines.append(f"beta = {self.beta}")
            lines.append(f"pH = {self.pH}")
            lines.append("")

        # Species definitions
        lines.append("species_definitions = {}")

        # Gas pressures
        for name, gas in self.gases.items():
            lines.append(f"species_definitions['{name}_g'] = {{'pressure': {gas.pressure}}}")

        # Add proton-electron if electrochemical
        if self.voltage is not None:
            lines.append("species_definitions['pe_g'] = {'pressure': 1.0}")

        # Site definition
        lines.append(f"species_definitions['{self.site_name}'] = {{'site_names': ['{self.site_type}'], 'total': {self.total_sites}}}")
        lines.append("")

        # Files
        lines.append(f"data_file = 'model.pkl'")
        lines.append(f"input_file = '{os.path.basename(energies_file)}'")
        lines.append("")

        # Thermodynamic modes
        lines.append(f"gas_thermo_mode = '{self.gas_thermo_mode}'")
        lines.append(f"adsorbate_thermo_mode = '{self.adsorbate_thermo_mode}'")
        if self.voltage is not None:
            lines.append(f"electrochemical_thermo_mode = '{self.electrochemical_thermo_mode}'")
        lines.append("")

        # Scaling constraints
        if self.scaling_constraints:
            lines.append("scaling_constraint_dict = {")
            for species, constraint in self.scaling_constraints.items():
                if isinstance(constraint, str):
                    lines.append(f"    '{species}': '{constraint}',")
                else:
                    lines.append(f"    '{species}': {constraint},")
            lines.append("}")
        else:
            # Auto-generate default constraints based on descriptor order
            lines.append("scaling_constraint_dict = {")
            # Get descriptor base names (without _s)
            desc_base_names = [d.replace('_s', '') for d in self.descriptor_names] if self.descriptor_names else []

            for i, name in enumerate(self.adsorbates.keys()):
                constraint = [0] * len(desc_base_names) + [None]
                if name in desc_base_names:
                    idx = desc_base_names.index(name)
                    constraint[idx] = '+'
                lines.append(f"    '{name}_s': {constraint},")

            for ts_name, ts in self.transition_states.items():
                mode = ts.scaling_mode
                if isinstance(mode, str):
                    lines.append(f"    '{ts_name}_s': '{mode}',")
                else:
                    lines.append(f"    '{ts_name}_s': {mode},")
            lines.append("}")
        lines.append("")

        # Interaction model
        if self.adsorbate_interaction_model:
            lines.append(f"adsorbate_interaction_model = '{self.adsorbate_interaction_model}'")
            lines.append("interaction_response_function = 'smooth_piecewise_linear'")
            lines.append(f"species_definitions['{self.site_name}']['interaction_response_parameters'] = {{'cutoff': 0.25, 'smoothing': 0.01}}")
            for sp, params in self.interaction_parameters.items():
                for param_type, values in params.items():
                    lines.append(f"species_definitions['{sp}_s'] = {{'{param_type}': {values}}}")
            lines.append("")

        # Atomic reservoir list (if set) - used in _generate_mkm_file
        if self.atomic_reservoir_list:
            lines.append(f"atomic_reservoir_list = {self.atomic_reservoir_list}")
            lines.append("")

        # Solver settings
        lines.append(f"decimal_precision = {self.decimal_precision}")
        lines.append(f"tolerance = {self.tolerance}")
        lines.append(f"max_rootfinding_iterations = {self.max_iterations}")
        lines.append(f"max_bisections = {self.max_bisections}")
        if self.use_numbers_solver is not None:
            lines.append(f"use_numbers_solver = {self.use_numbers_solver}")
        if self.max_damping_iterations is not None:
            lines.append(f"max_damping_iterations = {self.max_damping_iterations}")
        if self.max_initial_guesses is not None:
            lines.append(f"max_initial_guesses = {self.max_initial_guesses}")
        lines.append("")

        with open(filepath, 'w') as f:
            f.write('\n'.join(lines))

    # =========================================================================
    # Running the model
    # =========================================================================

    def run(self, output_dir: Optional[str] = None) -> 'CatMAPResult':
        """
        Run the microkinetic model.

        Args:
            output_dir: Directory for output files (uses work_dir if None)

        Returns:
            CatMAPResult object containing all results
        """
        from catmap.model import ReactionModel
        from catmap import analyze

        output_dir = output_dir or self.work_dir
        os.makedirs(output_dir, exist_ok=True)

        # Generate input files
        energies_file = os.path.join(output_dir, 'energies.txt')
        mkm_file = os.path.join(output_dir, 'model.mkm')

        self._generate_energies_file(energies_file)
        self._generate_mkm_file(mkm_file, energies_file)

        # Change to output directory and run
        original_dir = os.getcwd()
        try:
            os.chdir(output_dir)

            # Create and run model
            self.model = ReactionModel(setup_file='model.mkm')
            self.model.output_variables = self.output_variables
            self.model.run()

            # Create result object
            self.results = CatMAPResult(
                model=self.model,
                output_dir=output_dir,
                predictor=self
            )

        finally:
            os.chdir(original_dir)

        return self.results

    def run_single_point(self, surface: str, output_dir: Optional[str] = None) -> 'SinglePointResult':
        """
        Run microkinetic model for a single surface (single point calculation).

        This is useful when you want to calculate rates/coverages for a specific
        catalyst without scanning the entire descriptor space.

        Args:
            surface: Surface name (must have energies defined)
            output_dir: Directory for output files

        Returns:
            SinglePointResult object with rates and coverages
        """
        from catmap.model import ReactionModel

        output_dir = output_dir or self.work_dir
        os.makedirs(output_dir, exist_ok=True)

        # Generate input files
        energies_file = os.path.join(output_dir, 'energies.txt')
        mkm_file = os.path.join(output_dir, 'single_point.mkm')

        self._generate_energies_file(energies_file)
        self._generate_single_point_mkm(mkm_file, energies_file, surface)

        # Change to output directory and run
        original_dir = os.getcwd()
        try:
            os.chdir(output_dir)

            # Create and run model
            model = ReactionModel(setup_file='single_point.mkm')
            model.output_variables = ['coverage', 'rate', 'production_rate', 'free_energy']
            model.run()

            # Extract results at the single point
            result = SinglePointResult(
                surface=surface,
                model=model,
                output_dir=output_dir,
                predictor=self
            )

        finally:
            os.chdir(original_dir)

        return result

    def _generate_single_point_mkm(self, filepath: str, energies_file: str, surface: str):
        """Generate mkm file for single point calculation."""
        lines = []

        # Reaction expressions
        lines.append("rxn_expressions = [")
        for rxn in self.reactions:
            lines.append(f"    '{rxn}',")
        lines.append("]")
        lines.append("")

        # Only the specified surface
        lines.append(f"surface_names = ['{surface}']")
        lines.append("")

        # Get descriptors for this surface from adsorbate energies
        desc_values = []
        for desc in self.descriptor_names:
            desc_base = desc.replace('_s', '')
            if desc_base in self.adsorbates:
                energy = self.adsorbates[desc_base].formation_energies.get(surface)
                if energy is not None:
                    desc_values.append(energy)

        if not desc_values:
            # Auto-detect from first two adsorbates
            ads_list = list(self.adsorbates.keys())[:2]
            desc_values = []
            desc_names = []
            for ads in ads_list:
                energy = self.adsorbates[ads].formation_energies.get(surface)
                if energy is not None:
                    desc_values.append(energy)
                    desc_names.append(f"{ads}_s")
            self.descriptor_names = desc_names

        # Single point: resolution=1, range is just the descriptor value
        lines.append(f"descriptor_names = {self.descriptor_names}")
        ranges = [[v, v] for v in desc_values]
        lines.append(f"descriptor_ranges = {ranges}")
        lines.append("resolution = 1")
        lines.append("")

        # Temperature
        lines.append(f"temperature = {self.temperature}")
        lines.append("")

        # Electrochemical settings
        if self.voltage is not None:
            lines.append(f"voltage = {self.voltage}")
            lines.append(f"beta = {self.beta}")
            lines.append(f"pH = {self.pH}")
            lines.append("")

        # Species definitions
        lines.append("species_definitions = {}")
        for name, gas in self.gases.items():
            lines.append(f"species_definitions['{name}_g'] = {{'pressure': {gas.pressure}}}")
        if self.voltage is not None:
            lines.append("species_definitions['pe_g'] = {'pressure': 1.0}")
        lines.append(f"species_definitions['{self.site_name}'] = {{'site_names': ['{self.site_type}'], 'total': {self.total_sites}}}")
        lines.append("")

        # Ideal-gas parameters (symmetry, geometry, spin) for any gas whose
        # CatMAP built-in table would miss it. Must be set before load() so
        # ideal_gas thermo can find it.
        ideal_gas_entries = {}
        atoms_entries = []  # (species_key, symbols_list, positions_nested_list)
        for name, gas in self.gases.items():
            if gas.symmetry_number is None or gas.geometry is None:
                continue
            spin = 0 if gas.spin is None else gas.spin
            ideal_gas_entries[f"{name}_g"] = [gas.symmetry_number, gas.geometry, spin]
            if gas.atoms_geometry is not None:
                try:
                    symbols = list(gas.atoms_geometry.get_chemical_symbols())
                    positions = [
                        [float(p[0]), float(p[1]), float(p[2])]
                        for p in gas.atoms_geometry.get_positions()
                    ]
                    atoms_entries.append((f"{name}_g", symbols, positions))
                except Exception:
                    pass
        if ideal_gas_entries:
            lines.append(
                "ideal_gas_params = "
                + repr(ideal_gas_entries)
            )
            lines.append("")
        if atoms_entries:
            # Inline Atoms construction: the mkm file is exec'd as Python, so
            # CatMAP picks up this atoms_dict before the thermo init runs.
            lines.append("from ase import Atoms as _ASE_Atoms")
            lines.append("atoms_dict = {}")
            for key, symbols, positions in atoms_entries:
                lines.append(
                    f"atoms_dict['{key}'] = _ASE_Atoms(symbols={symbols!r}, "
                    f"positions={positions!r})"
                )
            lines.append("")

        # Files
        lines.append("data_file = 'single_point.pkl'")
        lines.append(f"input_file = '{os.path.basename(energies_file)}'")
        lines.append("")

        # Thermodynamic modes
        lines.append(f"gas_thermo_mode = '{self.gas_thermo_mode}'")
        lines.append(f"adsorbate_thermo_mode = '{self.adsorbate_thermo_mode}'")
        if self.voltage is not None:
            lines.append(f"electrochemical_thermo_mode = '{self.electrochemical_thermo_mode}'")
        lines.append("")

        # Scaling - for single surface, use direct values
        lines.append("scaling_constraint_dict = {")
        desc_base_names = [d.replace('_s', '') for d in self.descriptor_names]
        for name in self.adsorbates.keys():
            constraint = [0] * len(desc_base_names) + [None]
            if name in desc_base_names:
                idx = desc_base_names.index(name)
                constraint[idx] = '+'
            lines.append(f"    '{name}_s': {constraint},")
        for ts_name, ts in self.transition_states.items():
            mode = ts.scaling_mode
            if isinstance(mode, str):
                lines.append(f"    '{ts_name}_s': '{mode}',")
            else:
                lines.append(f"    '{ts_name}_s': {mode},")
        lines.append("}")
        lines.append("")

        # Atomic reservoir list (if set)
        if self.atomic_reservoir_list:
            lines.append(f"atomic_reservoir_list = {self.atomic_reservoir_list}")
            lines.append("")

        # Solver settings
        lines.append(f"decimal_precision = {self.decimal_precision}")
        lines.append(f"tolerance = {self.tolerance}")
        lines.append(f"max_rootfinding_iterations = {self.max_iterations}")
        lines.append(f"max_bisections = {self.max_bisections}")
        if self.use_numbers_solver is not None:
            lines.append(f"use_numbers_solver = {self.use_numbers_solver}")
        if self.max_damping_iterations is not None:
            lines.append(f"max_damping_iterations = {self.max_damping_iterations}")
        if self.max_initial_guesses is not None:
            lines.append(f"max_initial_guesses = {self.max_initial_guesses}")
        lines.append("")

        with open(filepath, 'w') as f:
            f.write('\n'.join(lines))


# =============================================================================
# Single Point Result Class
# =============================================================================

class SinglePointResult:
    """Container for single point calculation results."""

    def __init__(self, surface: str, model, output_dir: str, predictor):
        self.surface = surface
        self.model = model
        self.output_dir = output_dir
        self.predictor = predictor
        self.temperature = getattr(predictor, "temperature", None)
        self.voltage = getattr(predictor, "voltage", None)

        # Extract results
        self._extract_results()

    def __getstate__(self):
        # ReactionModel and predictor runtime objects are not stable pickle payloads.
        # Persist only the extracted single-point observables plus minimal metadata.
        return {
            "surface": self.surface,
            "output_dir": self.output_dir,
            "temperature": self.temperature,
            "voltage": self.voltage,
            "coverages": self.coverages,
            "rates": self.rates,
            "production_rates": self.production_rates,
        }

    def __setstate__(self, state):
        self.surface = state.get("surface")
        self.output_dir = state.get("output_dir")
        self.temperature = state.get("temperature")
        self.voltage = state.get("voltage")
        self.coverages = state.get("coverages", {})
        self.rates = state.get("rates", {})
        self.production_rates = state.get("production_rates", {})
        self.model = None
        self.predictor = None

    def _extract_results(self):
        """Extract coverages and rates from the model."""
        # Get coverage at the single point
        coverage_map = getattr(self.model, 'coverage_map', None)
        rate_map = getattr(self.model, 'rate_map', None)
        production_rate_map = getattr(self.model, 'production_rate_map', None)

        self.coverages = {}
        self.rates = {}
        self.production_rates = {}

        if coverage_map and len(coverage_map) > 0:
            # coverage_map is list of [descriptors, [coverages]]
            point = coverage_map[0]
            coverage_values = point[1] if len(point) > 1 else []

            # Map to species names
            adsorbate_names = self.model.adsorbate_names
            for i, name in enumerate(adsorbate_names):
                if i < len(coverage_values):
                    self.coverages[name] = float(coverage_values[i])

        if rate_map and len(rate_map) > 0:
            point = rate_map[0]
            rate_values = point[1] if len(point) > 1 else []

            # Map to reaction names
            rxn_expressions = self.model.rxn_expressions
            for i, rxn in enumerate(rxn_expressions):
                if i < len(rate_values):
                    self.rates[f"rxn_{i+1}"] = float(rate_values[i])

        if production_rate_map and len(production_rate_map) > 0:
            point = production_rate_map[0]
            prod_values = point[1] if len(point) > 1 else []

            # Map to gas species
            gas_names = self.model.gas_names
            for i, name in enumerate(gas_names):
                if i < len(prod_values):
                    self.production_rates[name] = float(prod_values[i])

    def get_coverage(self, species: Optional[str] = None) -> Union[Dict, float]:
        """
        Get surface coverage.

        Args:
            species: Species name (e.g., 'CO_s'). None returns all.

        Returns:
            Coverage value or dict of all coverages
        """
        if species:
            if species in self.coverages:
                return self.coverages[species]
            if species.endswith("_s"):
                return self.coverages.get(species[:-2], 0.0)
            return self.coverages.get(f"{species}_s", self.coverages.get(species, 0.0))
        return self.coverages

    def get_rate(self, reaction: Optional[int] = None) -> Union[Dict, float]:
        """
        Get reaction rate.

        Args:
            reaction: Reaction index (1-indexed). None returns all.

        Returns:
            Rate value or dict of all rates
        """
        if reaction:
            return self.rates.get(f"rxn_{reaction}", 0.0)
        return self.rates

    def get_production_rate(self, species: Optional[str] = None) -> Union[Dict, float]:
        """
        Get production rate.

        Args:
            species: Gas species name (e.g., 'CO2_g'). None returns all.

        Returns:
            Production rate or dict of all rates
        """
        if species:
            if species in self.production_rates:
                return self.production_rates[species]
            if species.endswith("_g"):
                return self.production_rates.get(species[:-2], 0.0)
            return self.production_rates.get(f"{species}_g", self.production_rates.get(species, 0.0))
        return self.production_rates

    def get_turnover_frequency(self) -> float:
        """Get the overall turnover frequency (max production rate)."""
        if self.production_rates:
            return max(abs(v) for v in self.production_rates.values())
        return 0.0

    def summary(self) -> str:
        """Generate a text summary of results."""
        lines = ["=" * 60]
        lines.append(f"Single Point Results: {self.surface}")
        lines.append("=" * 60)
        if self.temperature is not None:
            lines.append(f"Temperature: {self.temperature} K")
        if self.voltage is not None:
            lines.append(f"Voltage: {self.voltage} V vs SHE")

        lines.append("")
        lines.append("Coverages:")
        for species, coverage in self.coverages.items():
            lines.append(f"  {species}: {coverage:.6f}")

        lines.append("")
        lines.append("Reaction Rates (1/s):")
        for rxn, rate in self.rates.items():
            lines.append(f"  {rxn}: {rate:.6e}")

        lines.append("")
        lines.append("Production Rates (1/s):")
        for species, rate in self.production_rates.items():
            lines.append(f"  {species}: {rate:.6e}")

        lines.append("")
        lines.append(f"Turnover Frequency: {self.get_turnover_frequency():.6e} 1/s")
        lines.append("=" * 60)
        return "\n".join(lines)


# =============================================================================
# Result Class
# =============================================================================

class CatMAPResult:
    """Container for CatMAP results with analysis methods."""

    def __init__(self, model, output_dir: str, predictor: CatMAPPredictor):
        self.model = model
        self.output_dir = output_dir
        self.predictor = predictor

    def get_coverage_map(self):
        """Get coverage map data."""
        return getattr(self.model, 'coverage_map', None)

    def get_rate_map(self):
        """Get rate map data."""
        return getattr(self.model, 'rate_map', None)

    def get_production_rate_map(self):
        """Get production rate map data."""
        return getattr(self.model, 'production_rate_map', None)

    def get_descriptor_dict(self):
        """Get descriptor dictionary mapping surfaces to descriptor values."""
        return getattr(self.model, 'descriptor_dict', None)

    def plot_rate(self, save: str = 'rate.png', **kwargs):
        """Plot reaction rate volcano."""
        from catmap import analyze
        vm = analyze.VectorMap(self.model)
        vm.plot_variable = 'rate'
        vm.log_scale = kwargs.get('log_scale', True)
        vm.min = kwargs.get('min', 1e-25)
        vm.max = kwargs.get('max', 1e5)
        save_path = os.path.join(self.output_dir, save)
        return vm.plot(save=save_path)

    def plot_coverage(self, save: str = 'coverage.png', species: Optional[List[str]] = None, **kwargs):
        """
        Plot surface coverage.

        Args:
            save: Output filename
            species: List of species to include (e.g., ['CO_s', 'O_s']). None for all.
        """
        from catmap import analyze
        vm = analyze.VectorMap(self.model)
        vm.plot_variable = 'coverage'
        vm.log_scale = kwargs.get('log_scale', False)
        vm.min = kwargs.get('min', 0)
        vm.max = kwargs.get('max', 1)
        if species:
            vm.include_labels = species
        save_path = os.path.join(self.output_dir, save)
        return vm.plot(save=save_path)

    def plot_production_rate(self, save: str = 'production_rate.png', **kwargs):
        """Plot production rate."""
        from catmap import analyze
        vm = analyze.VectorMap(self.model)
        vm.plot_variable = 'production_rate'
        vm.log_scale = kwargs.get('log_scale', True)
        vm.min = kwargs.get('min', 1e-25)
        vm.max = kwargs.get('max', 1e5)
        vm.threshold = kwargs.get('threshold', 1e-30)
        save_path = os.path.join(self.output_dir, save)
        return vm.plot(save=save_path)

    def plot_scaling(self, save: str = 'scaling.png'):
        """Plot scaling relations."""
        from catmap import analyze
        sa = analyze.ScalingAnalysis(self.model)
        save_path = os.path.join(self.output_dir, save)
        return sa.plot(save=save_path)

    def plot_mechanism(self, surfaces: Optional[List[str]] = None,
                       mechanism: Optional[List[int]] = None,
                       save: str = 'mechanism.png', **kwargs):
        """
        Plot potential energy diagram for reaction mechanism.

        Args:
            surfaces: List of surfaces to plot (default: all)
            mechanism: List of reaction step indices (1-indexed, negative for reverse)
            save: Output filename
        """
        from catmap import analyze
        ma = analyze.MechanismAnalysis(self.model)

        # Configure
        ma.energy_type = kwargs.get('energy_type', 'free_energy')
        ma.include_labels = kwargs.get('include_labels', True)
        ma.pressure_correction = kwargs.get('pressure_correction', True)

        save_path = os.path.join(self.output_dir, save)
        return ma.plot(plot_variants=surfaces, save=save_path)

    def plot_rate_control(self, save: str = 'rate_control.png', **kwargs):
        """Plot degree of rate control."""
        from catmap import analyze

        # Need to recalculate with rate_control output
        if 'rate_control' not in self.model.output_variables:
            print("Warning: rate_control not in output_variables. Recalculating...")
            self.model.output_variables.append('rate_control')
            self.model.run()

        mm = analyze.MatrixMap(self.model)
        mm.plot_variable = 'rate_control'
        mm.log_scale = kwargs.get('log_scale', False)
        mm.min = kwargs.get('min', -2)
        mm.max = kwargs.get('max', 2)
        save_path = os.path.join(self.output_dir, save)
        return mm.plot(save=save_path)

    def plot_selectivity(self, save: str = 'selectivity.png', **kwargs):
        """Plot product selectivity."""
        from catmap import analyze

        if 'selectivity' not in self.model.output_variables:
            print("Warning: selectivity not in output_variables. Recalculating...")
            self.model.output_variables.append('selectivity')
            self.model.run()

        vm = analyze.VectorMap(self.model)
        vm.plot_variable = 'selectivity'
        vm.log_scale = kwargs.get('log_scale', False)
        vm.min = kwargs.get('min', 0)
        vm.max = kwargs.get('max', 1)
        save_path = os.path.join(self.output_dir, save)
        return vm.plot(save=save_path)

    def plot_all(self, prefix: str = ''):
        """Generate all standard plots."""
        plots = {}
        try:
            plots['rate'] = self.plot_rate(save=f'{prefix}rate.png')
        except Exception as e:
            print(f"Warning: Could not plot rate: {e}")
        try:
            plots['coverage'] = self.plot_coverage(save=f'{prefix}coverage.png')
        except Exception as e:
            print(f"Warning: Could not plot coverage: {e}")
        try:
            plots['production_rate'] = self.plot_production_rate(save=f'{prefix}production_rate.png')
        except Exception as e:
            print(f"Warning: Could not plot production_rate: {e}")
        try:
            plots['scaling'] = self.plot_scaling(save=f'{prefix}scaling.png')
        except Exception as e:
            print(f"Warning: Could not plot scaling: {e}")
        return plots

    def get_rate_at_surface(self, surface: str) -> Dict:
        """
        Get reaction rates at a specific surface.

        Args:
            surface: Surface name (e.g., 'Pt')

        Returns:
            Dict with rates for each reaction
        """
        desc_dict = self.get_descriptor_dict()
        if desc_dict and surface in desc_dict:
            descriptors = desc_dict[surface]
            # Find nearest point in rate_map
            rate_map = self.get_rate_map()
            if rate_map:
                # rate_map is list of [descriptors, rates]
                min_dist = float('inf')
                best_rates = None
                for point in rate_map:
                    desc, rates = point[0], point[1]
                    dist = sum((d1 - d2)**2 for d1, d2 in zip(desc, descriptors))
                    if dist < min_dist:
                        min_dist = dist
                        best_rates = rates
                return {'surface': surface, 'descriptors': descriptors, 'rates': best_rates}
        return None

    def summary(self) -> str:
        """Generate a text summary of results."""
        lines = ["=" * 60]
        lines.append("CatMAP Microkinetic Model Results")
        lines.append("=" * 60)
        lines.append(f"Temperature: {self.predictor.temperature} K")
        if self.predictor.voltage is not None:
            lines.append(f"Voltage: {self.predictor.voltage} V vs SHE")
            lines.append(f"pH: {self.predictor.pH}")
        lines.append(f"Surfaces: {self.predictor.surfaces}")
        lines.append(f"Reactions: {len(self.predictor.reactions)}")
        for i, rxn in enumerate(self.predictor.reactions):
            lines.append(f"  {i+1}. {rxn}")
        lines.append(f"Adsorbates: {list(self.predictor.adsorbates.keys())}")
        lines.append(f"Transition States: {list(self.predictor.transition_states.keys())}")
        lines.append(f"Output Directory: {self.output_dir}")
        lines.append("=" * 60)
        return "\n".join(lines)


# =============================================================================
# Convenience functions for common reactions
# =============================================================================

def create_co_oxidation_model(energies: Dict[str, Dict[str, float]],
                               temperature: float = 500,
                               co_pressure: float = 1.0,
                               o2_pressure: float = 0.33) -> CatMAPPredictor:
    """
    Create a CO oxidation microkinetic model.

    Args:
        energies: Dict with keys 'CO', 'O', 'O-CO', 'O-O' mapping to {surface: energy}
        temperature: Reaction temperature in K
        co_pressure: CO partial pressure
        o2_pressure: O2 partial pressure

    Returns:
        Configured CatMAPPredictor
    """
    predictor = CatMAPPredictor()

    # Add reactions
    predictor.add_reaction('*_s + CO_g -> CO*')
    predictor.add_reaction('2*_s + O2_g <-> O-O* + *_s -> 2O*')
    predictor.add_reaction('CO* + O* <-> O-CO* + * -> CO2_g + 2*')

    # Add gases
    predictor.add_gas('CO', formation_energy=0.0, pressure=co_pressure, frequencies=[2170])
    predictor.add_gas('O2', formation_energy=0.0, pressure=o2_pressure, frequencies=[1580])
    predictor.add_gas('CO2', formation_energy=-2.45, pressure=0, frequencies=[1333, 2349, 667, 667])

    # Add adsorbates and TS
    predictor.add_adsorbate('CO', energies.get('CO', {}))
    predictor.add_adsorbate('O', energies.get('O', {}))
    predictor.add_transition_state('O-CO', energies.get('O-CO', {}), scaling_mode='initial_state')
    predictor.add_transition_state('O-O', energies.get('O-O', {}), scaling_mode='final_state')

    predictor.set_conditions(temperature=temperature)

    return predictor


def create_her_model(energies: Dict[str, Dict[str, float]],
                     temperature: float = 300,
                     voltage: float = -0.1,
                     pH: float = 0) -> CatMAPPredictor:
    """
    Create a HER (Hydrogen Evolution Reaction) microkinetic model.

    Args:
        energies: Dict with keys 'H', 'H-H', 'pe-H' mapping to {surface: energy}
        temperature: Reaction temperature in K
        voltage: Electrode potential vs SHE
        pH: Solution pH

    Returns:
        Configured CatMAPPredictor
    """
    predictor = CatMAPPredictor()

    # Add reactions (Volmer-Heyrovsky-Tafel mechanism)
    predictor.add_reaction('pe_g + *_s <-> H*')                           # Volmer
    predictor.add_reaction('pe_g + H* <-> pe-H* <-> H2_g + *_s')          # Heyrovsky
    predictor.add_reaction('H* + H* <-> H-H* + *_s <-> H2_g + 2*_s')      # Tafel

    # Add gases
    predictor.add_gas('H2', formation_energy=0.0, pressure=0.01, frequencies=[4424])
    predictor.add_gas('H2O', formation_energy=0.0, pressure=0.03, frequencies=[3950, 3834, 1571])

    # Add adsorbates and TS
    predictor.add_adsorbate('H', energies.get('H', {}))
    predictor.add_transition_state('H-H', energies.get('H-H', {}), scaling_mode='initial_state')
    predictor.add_transition_state('pe-H', energies.get('pe-H', {}), scaling_mode='initial_state')

    predictor.set_conditions(temperature=temperature, voltage=voltage, pH=pH)
    predictor.set_thermo_modes(electrochemical_mode='hbond_electrochemical')

    return predictor


def create_orr_model(energies: Dict[str, Dict[str, float]],
                     temperature: float = 298.15,
                     voltage: float = 0.9,
                     pH: float = 0,
                     mechanism: str = 'associative') -> CatMAPPredictor:
    """
    Create an ORR (Oxygen Reduction Reaction) microkinetic model.

    Supports two mechanisms:
    - 'associative': O2 -> OOH -> O + H2O -> OH -> H2O (4e- pathway)
    - 'dissociative': O2 -> 2O -> 2OH -> 2H2O

    Args:
        energies: Dict with keys for intermediates mapping to {surface: energy}
            Required keys depend on mechanism:
            - associative: 'OH', 'O', 'OOH' (and optionally TS)
            - dissociative: 'OH', 'O', 'O-O'
        temperature: Reaction temperature in K
        voltage: Electrode potential vs SHE
        pH: Solution pH
        mechanism: 'associative' or 'dissociative'

    Returns:
        Configured CatMAPPredictor
    """
    predictor = CatMAPPredictor()

    if mechanism == 'associative':
        # 4-electron associative pathway (most common on Pt-group metals)
        # O2 + * + H+ + e- -> OOH*
        # OOH* + H+ + e- -> O* + H2O
        # O* + H+ + e- -> OH*
        # OH* + H+ + e- -> H2O + *
        predictor.add_reaction('O2_g + *_s + pe_g -> OOH*')
        predictor.add_reaction('OOH* + pe_g -> O* + H2O_g')
        predictor.add_reaction('O* + pe_g -> OH*')
        predictor.add_reaction('OH* + pe_g -> H2O_g + *_s')

        # Add gases
        predictor.add_gas('O2', formation_energy=0.0, pressure=0.21, frequencies=[1580])
        predictor.add_gas('H2O', formation_energy=0.0, pressure=0.03, frequencies=[3950, 3834, 1571])

        # Add adsorbates
        predictor.add_adsorbate('OOH', energies.get('OOH', {}))
        predictor.add_adsorbate('O', energies.get('O', {}))
        predictor.add_adsorbate('OH', energies.get('OH', {}))

        # Set descriptors (typically OH and O binding)
        predictor.descriptor_names = ['OH_s', 'O_s']

    else:  # dissociative
        # Dissociative pathway
        # O2 + 2* -> 2O*
        # O* + H+ + e- -> OH*
        # OH* + H+ + e- -> H2O + *
        predictor.add_reaction('O2_g + 2*_s <-> O-O* + *_s -> 2O*')
        predictor.add_reaction('O* + pe_g -> OH*')
        predictor.add_reaction('OH* + pe_g -> H2O_g + *_s')

        # Add gases
        predictor.add_gas('O2', formation_energy=0.0, pressure=0.21, frequencies=[1580])
        predictor.add_gas('H2O', formation_energy=0.0, pressure=0.03, frequencies=[3950, 3834, 1571])

        # Add adsorbates and TS
        predictor.add_adsorbate('O', energies.get('O', {}))
        predictor.add_adsorbate('OH', energies.get('OH', {}))
        if 'O-O' in energies:
            predictor.add_transition_state('O-O', energies.get('O-O', {}), scaling_mode='final_state')

        predictor.descriptor_names = ['OH_s', 'O_s']

    predictor.set_conditions(temperature=temperature, voltage=voltage, pH=pH)
    predictor.set_thermo_modes(
        gas_mode='ideal_gas',
        adsorbate_mode='harmonic_adsorbate',
        electrochemical_mode='hbond_electrochemical'
    )

    return predictor


def create_co2rr_model(energies: Dict[str, Dict[str, float]],
                       temperature: float = 298.15,
                       voltage: float = -0.8,
                       pH: float = 7,
                       product: str = 'CO') -> CatMAPPredictor:
    """
    Create a CO2RR (CO2 Reduction Reaction) microkinetic model.

    Supports different product pathways:
    - 'CO': CO2 -> COOH -> CO (2e- pathway)
    - 'HCOOH': CO2 -> OCHO -> HCOOH (2e- pathway)
    - 'CH4': Full 8e- pathway to methane (simplified)

    Args:
        energies: Dict with keys for intermediates mapping to {surface: energy}
            Required keys depend on product:
            - CO: 'COOH', 'CO'
            - HCOOH: 'OCHO', 'HCOOH' (or use 'OCHO' only)
            - CH4: 'COOH', 'CO', 'CHO', 'CH2O', 'CH3O', 'O', 'OH'
        temperature: Reaction temperature in K
        voltage: Electrode potential vs SHE
        pH: Solution pH
        product: 'CO', 'HCOOH', or 'CH4'

    Returns:
        Configured CatMAPPredictor
    """
    predictor = CatMAPPredictor()

    if product == 'CO':
        # 2e- pathway to CO
        # CO2 + * + H+ + e- -> COOH*
        # COOH* + H+ + e- -> CO* + H2O
        # CO* -> CO + *
        predictor.add_reaction('CO2_g + *_s + pe_g -> COOH*')
        predictor.add_reaction('COOH* + pe_g -> CO* + H2O_g')
        predictor.add_reaction('CO* -> CO_g + *_s')

        # Add gases
        predictor.add_gas('CO2', formation_energy=0.0, pressure=1.0, frequencies=[1333, 2349, 667, 667])
        predictor.add_gas('CO', formation_energy=0.0, pressure=0.0, frequencies=[2170])
        predictor.add_gas('H2O', formation_energy=0.0, pressure=0.03, frequencies=[3950, 3834, 1571])

        # Add adsorbates
        predictor.add_adsorbate('COOH', energies.get('COOH', {}))
        predictor.add_adsorbate('CO', energies.get('CO', {}))

        predictor.descriptor_names = ['CO_s', 'COOH_s']

    elif product == 'HCOOH':
        # 2e- pathway to formic acid/formate
        # CO2 + * + H+ + e- -> OCHO*
        # OCHO* + H+ + e- -> HCOOH + *
        predictor.add_reaction('CO2_g + *_s + pe_g -> OCHO*')
        predictor.add_reaction('OCHO* + pe_g -> HCOOH_g + *_s')

        # Add gases
        predictor.add_gas('CO2', formation_energy=0.0, pressure=1.0, frequencies=[1333, 2349, 667, 667])
        predictor.add_gas('HCOOH', formation_energy=0.0, pressure=0.0, frequencies=[3570, 2943, 1770, 1387, 1229, 1105, 625, 1033, 638])
        predictor.add_gas('H2O', formation_energy=0.0, pressure=0.03, frequencies=[3950, 3834, 1571])

        # Add adsorbates
        predictor.add_adsorbate('OCHO', energies.get('OCHO', {}))

        # For single adsorbate, use it as both descriptors or add H
        if 'H' in energies:
            predictor.add_adsorbate('H', energies.get('H', {}))
            predictor.descriptor_names = ['OCHO_s', 'H_s']
        else:
            predictor.descriptor_names = ['OCHO_s']

    elif product == 'CH4':
        # Simplified 8e- pathway to methane (via CO)
        # CO2 -> COOH -> CO -> CHO -> CH2O -> CH3O -> CH4 + O -> OH -> H2O
        predictor.add_reaction('CO2_g + *_s + pe_g -> COOH*')
        predictor.add_reaction('COOH* + pe_g -> CO* + H2O_g')
        predictor.add_reaction('CO* + pe_g -> CHO*')
        predictor.add_reaction('CHO* + pe_g -> CH2O*')
        predictor.add_reaction('CH2O* + pe_g -> CH3O*')
        predictor.add_reaction('CH3O* + pe_g -> CH4_g + O*')
        predictor.add_reaction('O* + pe_g -> OH*')
        predictor.add_reaction('OH* + pe_g -> H2O_g + *_s')

        # Add gases
        predictor.add_gas('CO2', formation_energy=0.0, pressure=1.0, frequencies=[1333, 2349, 667, 667])
        predictor.add_gas('CH4', formation_energy=0.0, pressure=0.0, frequencies=[2917, 1534, 1534, 3019, 3019, 3019, 1306, 1306, 1306])
        predictor.add_gas('H2O', formation_energy=0.0, pressure=0.03, frequencies=[3950, 3834, 1571])

        # Add adsorbates
        for ads in ['COOH', 'CO', 'CHO', 'CH2O', 'CH3O', 'O', 'OH']:
            if ads in energies:
                predictor.add_adsorbate(ads, energies.get(ads, {}))

        predictor.descriptor_names = ['CO_s', 'OH_s']

    predictor.set_conditions(temperature=temperature, voltage=voltage, pH=pH)
    predictor.set_thermo_modes(
        gas_mode='ideal_gas',
        adsorbate_mode='harmonic_adsorbate',
        electrochemical_mode='hbond_electrochemical'
    )

    return predictor


def create_nrr_model(energies: Dict[str, Dict[str, float]],
                     temperature: float = 298.15,
                     voltage: float = -0.5,
                     pH: float = 0,
                     mechanism: str = 'associative') -> CatMAPPredictor:
    """
    Create an NRR (Nitrogen Reduction Reaction) microkinetic model.

    Supports two mechanisms:
    - 'associative': Distal or alternating hydrogenation
    - 'dissociative': N2 dissociation followed by hydrogenation

    Args:
        energies: Dict with keys for intermediates mapping to {surface: energy}
            Typical keys: 'N2H', 'N2H2', 'N', 'NH', 'NH2', 'NH3'
        temperature: Reaction temperature in K
        voltage: Electrode potential vs SHE
        pH: Solution pH
        mechanism: 'associative' or 'dissociative'

    Returns:
        Configured CatMAPPredictor
    """
    predictor = CatMAPPredictor()

    if mechanism == 'dissociative':
        # Dissociative pathway
        # N2 + 2* -> 2N*
        # N* + H+ + e- -> NH*
        # NH* + H+ + e- -> NH2*
        # NH2* + H+ + e- -> NH3 + *
        predictor.add_reaction('N2_g + 2*_s -> 2N*')
        predictor.add_reaction('N* + pe_g -> NH*')
        predictor.add_reaction('NH* + pe_g -> NH2*')
        predictor.add_reaction('NH2* + pe_g -> NH3_g + *_s')

        # Add adsorbates
        for ads in ['N', 'NH', 'NH2']:
            if ads in energies:
                predictor.add_adsorbate(ads, energies.get(ads, {}))

        predictor.descriptor_names = ['N_s', 'NH_s']

    else:  # associative (distal)
        # Associative distal pathway
        # N2 + * -> N2*
        # N2* + H+ + e- -> N2H*
        # N2H* + H+ + e- -> N + NH3 (or further intermediates)
        predictor.add_reaction('N2_g + *_s -> N2*')
        predictor.add_reaction('N2* + pe_g -> N2H*')
        predictor.add_reaction('N2H* + pe_g -> N* + NH3_g')
        predictor.add_reaction('N* + pe_g -> NH*')
        predictor.add_reaction('NH* + pe_g -> NH2*')
        predictor.add_reaction('NH2* + pe_g -> NH3_g + *_s')

        # Add adsorbates
        for ads in ['N2', 'N2H', 'N', 'NH', 'NH2']:
            if ads in energies:
                predictor.add_adsorbate(ads, energies.get(ads, {}))

        predictor.descriptor_names = ['N_s', 'N2H_s'] if 'N2H' in energies else ['N_s', 'NH_s']

    # Add gases
    predictor.add_gas('N2', formation_energy=0.0, pressure=1.0, frequencies=[2358])
    predictor.add_gas('NH3', formation_energy=0.0, pressure=0.0, frequencies=[3337, 950, 3414, 3414, 1627, 1627])
    predictor.add_gas('H2O', formation_energy=0.0, pressure=0.03, frequencies=[3950, 3834, 1571])

    predictor.set_conditions(temperature=temperature, voltage=voltage, pH=pH)
    predictor.set_thermo_modes(
        gas_mode='ideal_gas',
        adsorbate_mode='harmonic_adsorbate',
        electrochemical_mode='simple_electrochemical'
    )

    return predictor


# =============================================================================
# Main entry point for testing
# =============================================================================

if __name__ == '__main__':
    print("CatMAP Predictor - Testing CO Oxidation Example")
    print("=" * 60)

    # Example: CO oxidation on various metal surfaces
    # Data from Falsig et al. and Angew. Chem. Int. Ed., 47, 4835 (2008)
    energies = {
        'CO': {'Pt': 1.7, 'Pd': 1.55, 'Cu': 2.58, 'Ag': 2.99, 'Au': 3.04, 'Rh': 1.34, 'Ru': 1.3, 'Ni': 1.63},
        'O': {'Pt': 1.62, 'Pd': 1.55, 'Cu': 1.07, 'Ag': 2.05, 'Au': 2.61, 'Rh': 0.55, 'Ru': -0.07, 'Ni': 0.35},
        'O-CO': {'Pt': 4.04, 'Pd': 4.2, 'Cu': 4.18, 'Ag': 5.05, 'Au': 5.74, 'Rh': 3.1, 'Ru': 2.53, 'Ni': 3.25},
        'O-O': {'Pt': 5.35, 'Pd': 5.34, 'Cu': 4.74, 'Ag': 5.98, 'Au': 7.22, 'Rh': 3.79, 'Ru': 3.34, 'Ni': 3.5},
    }

    predictor = create_co_oxidation_model(energies, temperature=500)

    # Set descriptors explicitly
    predictor.set_descriptors(
        names=['O_s', 'CO_s'],
        ranges=[[-1, 3], [-0.5, 4]],
        resolution=15
    )

    print(f"Surfaces: {predictor.surfaces}")
    print(f"Reactions: {predictor.reactions}")
    print(f"Adsorbates: {list(predictor.adsorbates.keys())}")
    print(f"Transition States: {list(predictor.transition_states.keys())}")

    # =========================================================================
    # Test 1: Full descriptor space scan (volcano plot)
    # =========================================================================
    print("\n" + "=" * 60)
    print("Test 1: Full Descriptor Space Scan")
    print("=" * 60)

    result = predictor.run(output_dir='output/catmap_test')
    print("\n" + result.summary())

    # Generate plots
    print("\nGenerating plots...")
    result.plot_all()
    print(f"Plots saved to: {result.output_dir}")

    # =========================================================================
    # Test 2: Single point calculation for Pt
    # =========================================================================
    print("\n" + "=" * 60)
    print("Test 2: Single Point Calculation for Pt")
    print("=" * 60)

    single_result = predictor.run_single_point('Pt', output_dir='output/catmap_single_pt')
    print("\n" + single_result.summary())

    # =========================================================================
    # Test 3: Compare multiple surfaces
    # =========================================================================
    print("\n" + "=" * 60)
    print("Test 3: Compare Multiple Surfaces")
    print("=" * 60)

    for surface in ['Pt', 'Pd', 'Cu', 'Rh']:
        try:
            sp_result = predictor.run_single_point(surface, output_dir=f'output/catmap_single_{surface.lower()}')
            tof = sp_result.get_turnover_frequency()
            print(f"{surface}: TOF = {tof:.2e} 1/s")
        except Exception as e:
            print(f"{surface}: Error - {e}")
