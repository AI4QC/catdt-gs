# pMuTT Predictor - Simplified Wrapper for Thermochemistry Calculations

## Introduction

`pMuTT_predictor.py` is a user-friendly wrapper around the pMuTT (Python Multiscale Thermochemistry Toolbox) library, designed to minimize user input while providing access to powerful thermochemistry calculations.

pMuTT is developed by Prof. Vlachos's group at the University of Delaware and specializes in:
- Converting DFT data to thermodynamic properties (NASA polynomials)
- Applying coverage-dependent models for heterogeneous catalysis
- Generating Chemkin/Cantera input files for microkinetic modeling
- Calculating reaction thermodynamics

## Installation

### 1. Install pMuTT and Dependencies

The pMuTT library has been installed in the `catdt` environment with all required dependencies:

```bash
# Already installed in catdt environment:
# - pmutt (1.4.17)
# - ase (3.23.0)
# - numpy, scipy, pandas, matplotlib
# - pymongo, dnspython, networkx
# - pygal, xlrd, more_itertools, openpyxl, PyYAML, pymatgen
```

### 2. Verify Installation

```bash
conda activate catdt
python -c "import pmutt; print('pMuTT version:', pmutt.__version__)"
```

## Quick Start

### Example 1: Generate NASA Polynomials from DFT Data

```python
from ase.build import molecule
from pmutt_predictor import pMuTTPredictor

# Create predictor
predictor = pMuTTPredictor(work_dir='./output')

# Add species from DFT calculations
predictor.add_species(
    name='H2O',
    atoms=molecule('H2O'),              # ASE Atoms object
    potentialenergy=-14.22,             # DFT energy (eV)
    vib_wavenumbers=[3825, 3710, 1582], # cm^-1
    phase='gas'
)

# Generate NASA polynomials
nasa_models = predictor.generate_nasa_polynomials()

# Calculate thermodynamic properties at 500 K
props = predictor.calculate_thermo_at_T('H2O', temperature=500.0)
print(f"H at 500K: {props['H']:.2f} kJ/mol")
print(f"S at 500K: {props['S']:.3f} kJ/mol/K")
print(f"G at 500K: {props['G']:.2f} kJ/mol")
```

### Example 2: One-Line NASA Generation

```python
from pmutt_predictor import quick_nasa_from_dft
from ase.build import molecule

# Generate NASA polynomial in one line
nasa = quick_nasa_from_dft(
    name='CH4',
    atoms=molecule('CH4'),
    potentialenergy=-24.06,
    vib_wavenumbers=[2917, 1534, 1534, 3019, 3019, 3019, 1306, 1306, 1306],
    phase='gas'
)

# Use immediately
H = nasa.get_H(T=500, units='kJ/mol')
```

### Example 3: Coverage Effects for Surface Species

```python
from ase.build import fcc111
from ase.build import add_adsorbate

predictor = pMuTTPredictor()

# Add surface species
slab = fcc111('Pt', size=(2, 2, 4), vacuum=10.0)
co_slab = slab.copy()
add_adsorbate(co_slab, 'C', height=2.0, position='ontop')
add_adsorbate(co_slab, 'O', height=3.1, position='ontop')

predictor.add_species(
    name='CO*',
    atoms=co_slab,
    potentialenergy=-380.50,
    vib_wavenumbers=[380, 350, 270, 2050],
    phase='surface'
)

predictor.add_species(
    name='O*',
    atoms=o_slab,
    potentialenergy=-367.20,
    vib_wavenumbers=[450, 430, 350],
    phase='surface'
)

# Add coverage effect: O* coverage affects CO* adsorption
predictor.add_coverage_effect(
    affected_species='CO*',
    coverage_species='O*',
    intervals=[0.0, 0.3, 0.7, 1.0],      # Coverage intervals (ML)
    slopes=[-8.0, -5.0, -2.0]            # Energy change (kcal/mol per ML)
)

# Generate models with coverage effects included
nasa_models = predictor.generate_nasa_polynomials()
```

### Example 4: Reaction Thermodynamics

```python
predictor = pMuTTPredictor()

# Add species (gas phase)
predictor.add_species(name='H2', atoms=molecule('H2'), ...)
predictor.add_species(name='O2', atoms=molecule('O2'), ...)
predictor.add_species(name='H2O', atoms=molecule('H2O'), ...)

# Generate NASA polynomials
predictor.generate_nasa_polynomials()

# Add reaction: H2 + 0.5 O2 -> H2O
predictor.add_reaction(
    reactants=['H2', 'O2'],
    products=['H2O'],
    reactants_stoich=[1.0, 0.5],
    products_stoich=[1.0]
)

# Calculate reaction properties at 500 K
rxn_df = predictor.calculate_reaction_properties(temperature=500.0)
print(rxn_df)
```

### Example 5: Export to Chemkin/Cantera

```python
predictor = pMuTTPredictor()

# ... add species and generate NASA polynomials ...

# Export to Chemkin format
predictor.export_to_chemkin('thermo.dat')

# Export to Cantera format (separates gas and surface phases)
predictor.export_to_cantera(
    filename='mechanism.yaml',
    gas_species=['CO', 'O2', 'CO2'],
    surface_species=['CO*', 'O*', '*']
)

# Export to Excel for analysis
predictor.export_to_excel(
    filename='thermo_data.xlsx',
    temperature_range=[300, 500, 700, 900]
)
```

## Complete Example: CO Oxidation on Pt(111)

See `test_pmutt_advanced.py` for a comprehensive example that demonstrates:
- Gas phase species (CO, O2, CO2)
- Surface species (CO*, O*, vacant sites)
- Coverage effects (O coverage affects CO adsorption)
- Surface reactions
- Export to Chemkin/Cantera formats

Run the example:
```bash
python test_pmutt_advanced.py
```

Output files will be generated in `./pmutt_advanced_output/`:
- `CO_oxidation_thermo.dat` - Chemkin format
- `CO_oxidation.yaml` - Cantera format
- `CO_oxidation_thermo.xlsx` - Excel spreadsheet

## Key Features

### 1. Simplified Interface
- **Minimal input required**: Only need DFT energy, vibrations, and structure
- **Automatic element detection**: From ASE Atoms objects
- **Automatic geometry detection**: Atom, linear, or nonlinear
- **Smart defaults**: Reasonable temperature ranges, symmetry numbers, etc.

### 2. Coverage Effects
- **Piecewise linear model**: Specify intervals and slopes
- **Multiple coverage effects**: Can have multiple species affecting each other
- **Automatic integration**: Coverage effects included in NASA polynomial generation

### 3. Flexible Input
- **ASE Atoms objects**: Standard structure format
- **VASP support**: Read from OUTCAR files
- **Gaussian support**: Read from log files
- **Manual input**: Specify all parameters directly

### 4. Multiple Output Formats
- **Chemkin**: Standard microkinetic modeling format
- **Cantera**: Modern Python-based reactor simulator
- **Excel**: Detailed thermodynamic tables
- **JSON**: For programmatic access

## Core Functionality Summary

### Species Management
| Method | Description |
|--------|-------------|
| `add_species()` | Add species with DFT data |
| `add_species_from_vasp()` | Read from VASP OUTCAR |
| `add_species_from_gaussian()` | Read from Gaussian log |
| `set_reference_species()` | Set DFT reference species |

### Coverage Effects
| Method | Description |
|--------|-------------|
| `add_coverage_effect()` | Add piecewise coverage model |

### Model Generation
| Method | Description |
|--------|-------------|
| `generate_statmech_models()` | Create statistical mechanics models |
| `generate_nasa_polynomials()` | Generate NASA polynomials |

### Reactions
| Method | Description |
|--------|-------------|
| `add_reaction()` | Add elementary reaction |
| `calculate_reaction_properties()` | Get ΔH, ΔS, ΔG for all reactions |

### Calculations
| Method | Description |
|--------|-------------|
| `calculate_thermo_at_T()` | Get Cp, H, S, G at temperature |

### Export
| Method | Description |
|--------|-------------|
| `export_to_chemkin()` | Export to Chemkin THERMO.DAT |
| `export_to_cantera()` | Export to Cantera YAML |
| `export_to_excel()` | Export to Excel spreadsheet |

### Utilities
| Method | Description |
|--------|-------------|
| `get_species_summary()` | Get DataFrame of all species |
| `summary()` | Print predictor status |

## Coverage Effect Models

The wrapper supports **piecewise linear coverage effects**, which is the "fast and practical" approach used in industry:

```
E_ads(θ) = E_0 + f(θ)

where f(θ) is piecewise linear:
  θ ∈ [0, 0.3]:    f(θ) = -8.0 * θ
  θ ∈ [0.3, 0.7]:  f(θ) = -8.0 * 0.3 + (-5.0) * (θ - 0.3)
  θ ∈ [0.7, 1.0]:  f(θ) = ... + (-2.0) * (θ - 0.7)
```

This model captures:
- **Strong repulsion at low coverage** (adsorbate-adsorbate interactions)
- **Moderate repulsion at medium coverage** (site blocking)
- **Weak repulsion at high coverage** (saturation effects)

## Thermodynamic Modes

The wrapper automatically selects appropriate thermodynamic models:

| Phase | Geometry | Translational | Rotational | Vibrational |
|-------|----------|---------------|------------|-------------|
| Gas   | Atom     | FreeTrans     | None       | None        |
| Gas   | Linear/Nonlinear | FreeTrans | RigidRotor | HarmonicVib |
| Surface | Any   | None          | None       | HarmonicVib |

## Temperature Ranges

Default NASA polynomial temperature ranges:
- **T_low**: 200 K
- **T_mid**: 1000 K
- **T_high**: 3500 K (gas), 2000 K (surface)

Customize with:
```python
nasa = predictor.generate_nasa_polynomials(
    T_low=298.15,
    T_mid=1000.0,
    T_high=2500.0
)
```

## Tips and Best Practices

### 1. DFT Energy References
- For gas phase species: Use consistent reference states (H2, H2O, N2, etc.)
- For surface species: Include full slab energy (not just adsorbate)
- Units: Always eV for DFT energies

### 2. Vibrational Frequencies
- Remove imaginary frequencies (use only real modes)
- For surface species: Include frustrated rotations/translations
- For transition states: Remove the reaction coordinate (imaginary mode)
- Units: Always cm^-1

### 3. Coverage Effects
- Start with simple linear models before using piecewise
- Use DFT calculations at different coverages to parameterize slopes
- Typical values: -5 to -10 kcal/mol per ML for repulsive interactions

### 4. Export Formats
- **Chemkin**: Best for legacy microkinetic codes (Chemkin, DETCHEM)
- **Cantera**: Best for Python-based reactor simulations
- **Excel**: Best for manual analysis and plotting

## Differences from CatMAP

While both pMuTT and CatMAP can be used for heterogeneous catalysis modeling, they have different focuses:

| Feature | pMuTT | CatMAP |
|---------|-------|--------|
| **Primary Use** | Generate NASA polynomials and thermo data | Mean-field microkinetic modeling |
| **Coverage Effects** | Semi-empirical linear models | Linear, piecewise, or cluster expansion |
| **Outputs** | NASA polynomials, Chemkin/Cantera files | Volcano plots, coverage maps, TOF |
| **Kinetics** | No kinetic solver (exports to Cantera/Chemkin) | Built-in mean-field solver |
| **Typical Workflow** | DFT → NASA → Cantera → reactor simulation | DFT → scaling relations → volcano plot |

**Recommended workflow**:
1. Use **pMuTT** to generate thermodynamic data from DFT
2. Use **CatMAP** for descriptor-based screening and volcano plots
3. Use **Cantera** (with pMuTT-generated data) for detailed reactor simulations

## Troubleshooting

### Issue: "Reference.__init__() missing required positional arguments"
**Solution**: References are optional. For basic usage, skip `set_reference_species()`.

### Issue: "AttributeError: ... object has no attribute ..."
**Solution**: Check that you've called `generate_nasa_polynomials()` before trying to export.

### Issue: Incorrect parameter order in get_Cp, get_H, etc.
**Solution**: Always use keyword arguments:
```python
# Correct
nasa.get_Cp(T=500, units='J/mol/K')

# Incorrect
nasa.get_Cp('J/mol/K', 500)  # Wrong order!
```

### Issue: Cantera export fails
**Solution**: Ensure phase is correctly specified ('gas' or 'surface') for all species.

## References

1. **pMuTT Documentation**: https://vlachosgroup.github.io/pMuTT/
2. **pMuTT Paper**: J. Lym et al., *Computer Physics Communications* (2019), DOI: 10.1016/j.cpc.2019.106864
3. **Cantera**: https://cantera.org/
4. **ASE**: https://wiki.fysik.dtu.dk/ase/

## Citation

If you use pMuTT in your research, please cite:

```
J. Lym, G.R. Wittreich and D.G. Vlachos,
A Python Multiscale Thermochemistry Toolbox (pMuTT) for thermochemical and kinetic parameter estimation,
Computer Physics Communications (2019) 106864,
https://doi.org/10.1016/j.cpc.2019.106864
```

## License

This wrapper script is part of CatDT and is licensed under the Apache-2.0 License.

pMuTT itself is licensed under the MIT License. Please refer to the pMuTT repository for details.

## Contact

For issues related to the wrapper, please check the CatDT repository.
For issues related to pMuTT itself, please visit: https://github.com/VlachosGroup/pMuTT/issues
