# UniMech — Reaction Mechanism Search

UniMech is CatDT's multi-pathway reaction mechanism discovery engine. It systematically explores reaction networks to find competing pathways from an initial state to a target product.

## Architecture

```
User Input (initial_state, target_state)
    ↓
CandidateGeneratorRouter
    └─ Generic RDKit Operations (element-agnostic)
    ↓
MechanismSearchEngine (beam search + sibling pruning)
    ↓
FreeEnergyRouter (UMA thermal / CHE electrochemical)
    ↓
Ranked pathways → Agent4/5 handoff or standalone NEB
```

## Three Exploration Modes

### 1. Agent-Guided (fast)
Uses LLM to recommend 3+ competing pathways, evaluates each with UMA.
Best for: well-known reactions, quick screening.

### 2. Systematic CRN (thorough)
Three-phase pipeline:
- Phase 1: BFS graph enumeration via RDKit bond operations (no energy calls, 30-50 raw paths)
- Phase 2: LLM chemical plausibility filter (select top 8-10)
- Phase 3: UMA energy evaluation on filtered paths

### 3. Combined (both_sequential / both_parallel)
Runs agent-guided first, then systematic for uncovered branches, or both in parallel.

## Candidate Generation — RDKit Operations

Element-agnostic bond operations that enumerate all reachable states:

| Operation | Description | Example |
|-----------|-------------|---------|
| Bond dissociation | Break a bond in adsorbate | *CH3OH → *CH3O + *H |
| Hydrogenation | Add H to adsorbate | *CO + *H → *CHO |
| Surface association | Combine two adsorbates | *CH3 + *O → *CH3O |
| Eley-Rideal | Adsorbate + gas molecule | *O + CO(g) → *CO2 |
| Associative desorption | Two adsorbates → gas product | 2*H → H2(g) |
| PCET | Proton-coupled electron transfer (electro) | *O + H⁺+e⁻ → *OH |
| Isomerization | Internal rearrangement | *CHOH → *CH2O |

## Key Classes & APIs

### MechanismSearchEngine
```python
from core.pathway.search_engine import MechanismSearchEngine

engine = MechanismSearchEngine(generator, energy_router, policy)
engine.initialize_root(species_label, elements_dict)
result = engine.search_until_stop(target_label, surface_path, backend)
# result: {"pathways": [...], "stats": {...}}
```

### CandidateGeneratorRouter
```python
from core.pathway.candidate_generators import CandidateGeneratorRouter

router = CandidateGeneratorRouter(enable_generic=True)
candidates = router.generate_candidates(mode="systematic", ...)
```

### SystematicCRNExplorer
```python
from core.pathway.candidate_generators import SystematicCRNExplorer

explorer = SystematicCRNExplorer(max_carbon=6, max_heavy_atoms=12)
children = explorer.expand_from_state(current_state, budget, depth)
```

### FreeEnergyRouter
```python
from core.pathway.free_energy_router import FreeEnergyRouter, ThermalStateEnergyBackend

thermal = ThermalStateEnergyBackend(tools=tools_instance)
router = FreeEnergyRouter(thermal_backend=thermal)
energy = router.evaluate(state, backend="thermal_uma")
```

### PruningPolicy
```python
from core.pathway.pruning_policy import PruningPolicy

policy = PruningPolicy(
    max_depth=8,
    beam_width=8,
    max_state_evaluations=40,
    delta_keep=0.10,   # eV: siblings within this gap → keep both
    delta_prune=0.30,  # eV: siblings beyond this gap → prune higher
)
```

## Via CatDTTools (Tool Facade — Recommended)

```python
from camel_agents.tools import CatDTTools

tools = CatDTTools(output_base_dir="output/mechanism")

# Step 1: Initialize context
ctx = tools.initialize_mechanism_context_tool(
    initial_state="*CO",
    target_state="CH4(g)",
    surface_id="Cu111",
    surface_structure_path="slab.vasp",
    clean_slab_path="clean_slab.vasp",
    environment={"T": 523, "P": 1e6},
)

# Step 2: Generate candidate steps
candidates = tools.generate_candidate_steps_tool(
    use_generic_ops=True, electro=False
)

# Step 3: Run search
result = tools.run_mechanism_search_tool(
    max_depth=8, beam_width=8,
    delta_keep=0.10, delta_prune=0.30,
    max_state_evaluations=40,
    energy_backend="thermal_uma",
    exploration_mode="both_sequential",
)

# Step 4: Export for NEB
shortlist = tools.extract_pathway_shortlist_tool(
    search_result=result,
    output_base_dir="output/mechanism",
    run_id="search_001",
)
```

## Data Schemas

### CatalyticStateRecord
```python
# from camel_agents.mechanism_schemas
CatalyticStateRecord(
    state_id="Cu111_hollow_CO*",
    species_label="*CO",
    phase=Phase.ADSORBED,  # ADSORBED | GAS | SOLVATED | SURFACE
    elements={"C": 1, "O": 1},
    surface_id="Cu111",
    site_type="hollow",     # hollow | bridge | atop
)
```

### ElementaryStepCandidate
```python
ElementaryStepCandidate(
    step_id="step_001",
    reactant_state_id="...",
    product_state_id="...",
    operation_type=OpType.DISSOCIATION,
    # ADSORPTION | DESORPTION | DISSOCIATION | ASSOCIATION |
    # HYDROGENATION | DEHYDROGENATION | PCET | COUPLING |
    # REARRANGEMENT | ELEY_RIDEAL | OTHER
    bond_changes="break C-H",
    estimated_barrier_eV=0.85,
    source_backend="generic_ops",  # generic_ops | llm
    confidence=0.8,
)
```

### CandidatePathway
```python
CandidatePathway(
    pathway_id="path_001",
    states=[...],          # ordered CatalyticStateRecord list
    steps=[...],           # ElementaryStepCandidate list
    total_free_energy_change_eV=-1.2,
    max_step_energy_eV=0.85,
    is_retained=True,
)
```

## Configuration (via CatDTConfig)

```python
CatDTConfig(
    enable_mechanism_search=True,
    mechanism_exploration_mode="agent_guided",  # or systematic, both_sequential, both_parallel
    mechanism_max_depth=8,
    mechanism_beam_width=8,
    mechanism_delta_keep=0.10,
    mechanism_delta_prune=0.30,
    mechanism_max_state_evaluations=40,
)
```

## Key Files

| File | Class | Purpose |
|------|-------|---------|
| `core/pathway/search_engine.py` | `MechanismSearchEngine` | Beam search algorithm |
| `core/pathway/candidate_generators.py` | `SystematicCRNExplorer`, `AgentGuidedPathwayGenerator` | Candidate generation |
| `core/pathway/free_energy_router.py` | `FreeEnergyRouter` | Energy evaluation routing |
| `core/pathway/pruning_policy.py` | `PruningPolicy` | Sibling comparison + beam selection |
| `core/pathway/pathway_bundle.py` | — | Export to Agent4/5 format |
| `camel_agents/mechanism_search.py` | `MechanismStageRunner` | High-level orchestrator |
| `camel_agents/mechanism_schemas.py` | Data models | State, step, pathway schemas |
| `camel_agents/tooling/mechanism.py` | `MechanismToolsMixin` | Tool implementations |
