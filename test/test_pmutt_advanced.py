#!/usr/bin/env python
"""
Advanced test script for pMuTT predictor demonstrating coverage effects
and surface species thermodynamics.

This script demonstrates:
1. Surface species thermodynamics
2. Coverage-dependent adsorption energies
3. Surface reactions
4. Export to Chemkin/Cantera formats
"""
import sys
sys.path.insert(0, "core/kmc")

from ase.build import molecule, fcc111
from ase.build import add_adsorbate
from pmutt_predictor import pMuTTPredictor

def test_surface_reactions_with_coverage():
    """
    Test CO oxidation on Pt(111) surface with coverage effects.

    Demonstrates:
    - Gas phase species (CO, O2, CO2)
    - Surface species (CO*, O*, vacant sites)
    - Coverage effects (O coverage affects CO adsorption)
    - Surface reactions
    """
    print("=" * 70)
    print("Advanced pMuTT Test: CO Oxidation with Coverage Effects")
    print("=" * 70)

    predictor = pMuTTPredictor(work_dir='output/pmutt', verbose=True)

    # =========================================================================
    # 1. Add gas phase species
    # =========================================================================
    print("\n" + "=" * 70)
    print("Step 1: Adding gas phase species")
    print("=" * 70)

    # CO gas
    predictor.add_species(
        name='CO',
        atoms=molecule('CO'),
        potentialenergy=-14.79,  # eV (example DFT energy)
        vib_wavenumbers=[2170],
        symmetrynumber=1,
        phase='gas'
    )

    # O2 gas
    predictor.add_species(
        name='O2',
        atoms=molecule('O2'),
        potentialenergy=-9.86,
        vib_wavenumbers=[1580],
        symmetrynumber=2,
        spin=1,
        phase='gas'
    )

    # CO2 gas
    predictor.add_species(
        name='CO2',
        atoms=molecule('CO2'),
        potentialenergy=-22.99,
        vib_wavenumbers=[1333, 2349, 667, 667],
        symmetrynumber=2,
        phase='gas'
    )

    # =========================================================================
    # 2. Add surface species
    # =========================================================================
    print("\n" + "=" * 70)
    print("Step 2: Adding surface species (adsorbates)")
    print("=" * 70)

    # Create Pt(111) surface for structure
    slab = fcc111('Pt', size=(2, 2, 4), vacuum=10.0)

    # CO* (adsorbed CO)
    co_slab = slab.copy()
    add_adsorbate(co_slab, 'C', height=2.0, position='ontop')
    add_adsorbate(co_slab, 'O', height=3.1, position='ontop')
    predictor.add_species(
        name='CO*',
        atoms=co_slab,
        potentialenergy=-380.50,  # eV (slab + CO)
        vib_wavenumbers=[380, 350, 270, 2050],  # Surface vibrations + CO stretch
        symmetrynumber=1,
        phase='surface',
        elements={'C': 1, 'O': 1}
    )

    # O* (adsorbed O)
    o_slab = slab.copy()
    add_adsorbate(o_slab, 'O', height=1.2, position='fcc')
    predictor.add_species(
        name='O*',
        atoms=o_slab,
        potentialenergy=-367.20,  # eV (slab + O)
        vib_wavenumbers=[450, 430, 350],
        symmetrynumber=1,
        phase='surface',
        elements={'O': 1}
    )

    # Empty site *
    predictor.add_species(
        name='*',
        atoms=slab,
        potentialenergy=-365.50,  # eV (clean slab)
        vib_wavenumbers=[],
        symmetrynumber=1,
        phase='surface',
        elements={}
    )

    # =========================================================================
    # 3. Add coverage effects
    # =========================================================================
    print("\n" + "=" * 70)
    print("Step 3: Adding coverage effects")
    print("=" * 70)
    print("O* coverage destabilizes CO* adsorption (piecewise linear model)")

    # O coverage affects CO adsorption
    # At low coverage: strong effect (-8 kcal/mol)
    # At medium coverage: moderate effect (-5 kcal/mol)
    # At high coverage: weak effect (-2 kcal/mol)
    predictor.add_coverage_effect(
        affected_species='CO*',
        coverage_species='O*',
        intervals=[0.0, 0.3, 0.7, 1.0],
        slopes=[-8.0, -5.0, -2.0]  # kcal/mol per ML
    )

    # Self-interaction: O coverage also affects O adsorption
    predictor.add_coverage_effect(
        affected_species='O*',
        coverage_species='O*',
        intervals=[0.0, 0.5, 1.0],
        slopes=[-3.0, -1.0]  # kcal/mol per ML
    )

    # =========================================================================
    # 4. Generate thermodynamic models
    # =========================================================================
    print("\n" + "=" * 70)
    print("Step 4: Generating NASA polynomials")
    print("=" * 70)

    nasa_models = predictor.generate_nasa_polynomials(
        T_low=200.0,
        T_mid=1000.0,
        T_high=2000.0  # Lower T_high for surface species
    )

    # =========================================================================
    # 5. Calculate thermodynamics at reaction conditions
    # =========================================================================
    print("\n" + "=" * 70)
    print("Step 5: Calculating thermodynamics at 500 K")
    print("=" * 70)

    T = 500.0
    for species_name in ['CO', 'O2', 'CO2', 'CO*', 'O*']:
        props = predictor.calculate_thermo_at_T(species_name, T, units='kJ/mol')
        print(f"\n{species_name:5s}: Cp={props['Cp']:6.2f} kJ/mol/K, "
              f"H={props['H']:8.2f} kJ/mol, "
              f"S={props['S']:6.3f} kJ/mol/K, "
              f"G={props['G']:8.2f} kJ/mol")

    # =========================================================================
    # 6. Add surface reactions
    # =========================================================================
    print("\n" + "=" * 70)
    print("Step 6: Adding surface reactions")
    print("=" * 70)

    # Reaction 1: CO adsorption
    # CO(g) + * -> CO*
    predictor.add_reaction(
        reactants=['CO', '*'],
        products=['CO*'],
        reactants_stoich=[1.0, 1.0],
        products_stoich=[1.0]
    )

    # Reaction 2: O2 dissociative adsorption
    # O2(g) + 2* -> 2O*
    predictor.add_reaction(
        reactants=['O2', '*'],
        products=['O*'],
        reactants_stoich=[1.0, 2.0],
        products_stoich=[2.0]
    )

    # Reaction 3: CO oxidation
    # CO* + O* -> CO2(g) + 2*
    predictor.add_reaction(
        reactants=['CO*', 'O*'],
        products=['CO2', '*'],
        reactants_stoich=[1.0, 1.0],
        products_stoich=[1.0, 2.0]
    )

    # =========================================================================
    # 7. Calculate reaction properties
    # =========================================================================
    print("\n" + "=" * 70)
    print("Step 7: Calculating reaction thermodynamics at 500 K")
    print("=" * 70)

    rxn_df = predictor.calculate_reaction_properties(temperature=500.0, units='kJ/mol')
    print("\nReaction properties:")
    print(rxn_df.to_string(index=False))

    # =========================================================================
    # 8. Export to different formats
    # =========================================================================
    print("\n" + "=" * 70)
    print("Step 8: Exporting to Chemkin and Cantera formats")
    print("=" * 70)

    # Export to Chemkin
    chemkin_file = predictor.export_to_chemkin('CO_oxidation_thermo.dat')
    print(f"\nChemkin file created: {chemkin_file}")

    # Export to Cantera (separate gas and surface phases)
    cantera_file = predictor.export_to_cantera(
        filename='CO_oxidation.yaml',
        gas_species=['CO', 'O2', 'CO2'],
        surface_species=['CO*', 'O*', '*']
    )
    print(f"Cantera file created: {cantera_file}")

    # Export to Excel for detailed analysis
    excel_file = predictor.export_to_excel(
        filename='CO_oxidation_thermo.xlsx',
        temperature_range=[300, 400, 500, 600, 700, 800]
    )
    print(f"Excel file created: {excel_file}")

    # =========================================================================
    # 9. Summary
    # =========================================================================
    print("\n" + predictor.summary())

    # Species summary
    print("\nSpecies Summary:")
    print(predictor.get_species_summary().to_string(index=False))

    print("\n" + "=" * 70)
    print("Advanced test completed successfully!")
    print("=" * 70)

    return predictor


def test_quick_nasa_generation():
    """
    Demonstrate quick NASA polynomial generation using convenience function.
    """
    print("\n\n" + "=" * 70)
    print("Quick NASA Generation Test")
    print("=" * 70)

    from pmutt_predictor import quick_nasa_from_dft

    # Generate NASA polynomial for CH4 in one line
    ch4_nasa = quick_nasa_from_dft(
        name='CH4',
        atoms=molecule('CH4'),
        potentialenergy=-24.06,
        vib_wavenumbers=[2917, 1534, 1534, 3019, 3019, 3019, 1306, 1306, 1306],
        phase='gas',
        symmetrynumber=12
    )

    # Calculate properties
    T = 500.0
    print(f"\nCH4 thermodynamics at {T} K:")
    print(f"  Cp = {ch4_nasa.get_Cp(T=T, units='J/mol/K'):.2f} J/mol/K")
    print(f"  H  = {ch4_nasa.get_H(T=T, units='kJ/mol'):.2f} kJ/mol")
    print(f"  S  = {ch4_nasa.get_S(T=T, units='J/mol/K'):.2f} J/mol/K")
    print(f"  G  = {ch4_nasa.get_G(T=T, units='kJ/mol'):.2f} kJ/mol")

    print("\n" + "=" * 70)


if __name__ == '__main__':
    # Run advanced test with coverage effects
    predictor = test_surface_reactions_with_coverage()

    # Run quick NASA generation test
    test_quick_nasa_generation()

    print("\n\n" + "=" * 70)
    print("ALL TESTS COMPLETED SUCCESSFULLY!")
    print("=" * 70)
    print("\nOutput files have been generated in:")
    print("  ./pmutt_advanced_output/")
    print("\nYou can now:")
    print("  1. Use the Chemkin file for microkinetic modeling")
    print("  2. Use the Cantera file for reactor simulations")
    print("  3. Analyze the Excel file for detailed thermodynamic data")
    print("=" * 70)
