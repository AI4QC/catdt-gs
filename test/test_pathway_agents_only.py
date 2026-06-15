from ase import Atom, Atoms


def test_align_added_atoms_to_anchor_preserves_offsets():
    from camel_agents import pathway_utils as pao

    anchor = [0.0, 0.0, 5.0]
    positions = [[1.0, 2.0, 6.0], [2.5, 4.0, 7.5]]
    aligned = pao.align_added_atoms_to_anchor(positions, anchor)

    assert aligned[0][0] == anchor[0]
    assert aligned[0][1] == anchor[1]
    assert aligned[0][2] == positions[0][2]

    dx = aligned[1][0] - aligned[0][0]
    dy = aligned[1][1] - aligned[0][1]
    dz = aligned[1][2] - aligned[0][2]
    assert dx == positions[1][0] - positions[0][0]
    assert dy == positions[1][1] - positions[0][1]
    assert dz == positions[1][2] - positions[0][2]


def test_select_anchor_position_prefers_given_index():
    from camel_agents import pathway_utils as pao

    slab = Atoms([Atom("Pt", (0, 0, 0)), Atom("Pt", (1, 0, 2.0))])
    anchor = pao.select_anchor_position(slab, surface_indices=[0, 1], anchor_index=1)
    assert anchor == [1.0, 0.0, 2.0]


def test_apply_step_modification_maps_adsorbate_indices():
    from camel_agents import pathway_utils as pao

    slab = Atoms([Atom("Pt", (0, 0, 0)), Atom("Pt", (1, 0, 0)), Atom("H", (0, 0, 1))])
    slab.info["surface_indices"] = [0, 1]
    slab.info["adsorbate_indices"] = [2]

    step = pao.PathwayStepSpec(
        step_name="*H -> *",
        reactant_formula="*H",
        product_formula="*",
        reaction_type="desorption",
        atoms_to_add=[],
        atoms_to_remove=[0],
        confidence=0.9,
    )

    product, adsorbate_indices = pao.apply_step_modification(
        slab,
        step,
        use_adsorbate_local_indices=True,
    )

    assert len(product) == 2
    assert adsorbate_indices == []


def test_apply_step_modification_ignores_surface_removal_when_no_adsorbate():
    from camel_agents import pathway_utils as pao

    slab = Atoms([Atom("Pt", (0, 0, 0)), Atom("Pt", (1, 0, 0))])
    slab.info["surface_indices"] = [0, 1]
    slab.info["adsorbate_indices"] = []

    step = pao.PathwayStepSpec(
        step_name="* -> *",
        reactant_formula="*",
        product_formula="*",
        reaction_type="noop",
        atoms_to_add=[],
        atoms_to_remove=[0],
        confidence=0.5,
    )

    product, adsorbate_indices = pao.apply_step_modification(
        slab,
        step,
        use_adsorbate_local_indices=True,
    )

    assert len(product) == 2
    assert adsorbate_indices == []


def test_default_renderer_is_not_matplotlib():
    from camel_agents import pathway_utils as pao

    assert pao.DEFAULT_RENDERER != "matplotlib"


def test_default_descriptions_include_stoichiometry_hints():
    from camel_agents import pathway_utils as pao

    her_desc = pao.get_default_her_description()
    nrr_desc = pao.get_default_nrr_description()

    assert "second hydrogen" in her_desc.lower() or "two hydrogen" in her_desc.lower()
    assert "two nh3" in nrr_desc.lower()


def test_enforce_minimum_height():
    from camel_agents import pathway_utils as pao

    positions = [[0.0, 0.0, 0.5], [1.0, 0.0, 2.0]]
    adjusted = pao.enforce_minimum_height(positions, surface_max_z=1.5, clearance=1.0)

    assert adjusted[0][2] >= 2.5
    assert adjusted[1][2] >= 2.5


def test_build_energy_diagram_inputs():
    from camel_agents import pathway_utils as pao

    intermediates = ["*", "*H", "*H2", "*"]
    energies = [0.0, -0.1, -0.2, 0.0]
    barriers = [0.2, 0.3, 0.1]

    labels, energy_list, barrier_map = pao.build_energy_diagram_inputs(
        intermediates, energies, barriers
    )

    assert labels == intermediates
    assert energy_list == energies
    assert barrier_map["* -> *H"] == 0.2
    assert barrier_map["*H -> *H2"] == 0.3
    assert barrier_map["*H2 -> *"] == 0.1
