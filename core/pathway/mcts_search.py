"""Monte Carlo Tree Search for mechanism pathway discovery.

Unlike the layer-wise beam/pruning search in ``search_engine.py``, MCTS
uses UCB1 to balance exploiting low-energy partial paths against exploring
branches that look bad early but may lead to lower-energy states later
("hard up front, easy after"). Intermediate energies are looked up
through a shared ``EnergyCache`` so the same species is never re-evaluated.

Reward convention (eV scale, shared with the rest of UniMech): all node
energies are corrected adsorption free energies ΔG_ads (see
``core.pathway.free_energy_router``), so rewards are O(1 eV) and the UCB1
exploration constant ``c_uct`` is meaningful. The reward of a simulation
is ``-max(ΔG_step)`` along the realised path when ``step_dg_fn`` is given
(stoichiometry-corrected step free energies — the preferred mode),
otherwise ``-max(ΔG_ads)`` over the path's evaluable nodes. Lower energy
= better = higher reward. Species that cannot be evaluated (None) simply
do not contribute to the bottleneck — they are unevaluable, never scored
as 0 from some other scale.
"""

from __future__ import annotations

import logging
import math
import random
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Set, Tuple

from core.pathway.candidate_generators import (
    CandidateGeneratorRouter,
    _check_rdkit,
    parse_species_elements,
)
from core.pathway.energy_cache import EnergyCache

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Generic valence-based plausibility check (MCTS-only, does NOT affect pruning)
# ---------------------------------------------------------------------------
# Applied only to candidates fetched by ``MCTSSearchEngine``. Pruning's own
# enumerate_all_paths_bfs path is unaffected and continues to see every
# candidate the underlying generator emits.


def _element_default_valence(symbol: str) -> Optional[int]:
    """Look up the default valence of ``symbol`` via RDKit's periodic table.

    Element-agnostic: covers every element RDKit knows about. Atoms without
    a defined default valence (e.g. most transition metals) return ``None``
    so callers can treat them permissively.
    """
    if not _check_rdkit():
        return None
    from rdkit import Chem
    try:
        pt = Chem.GetPeriodicTable()
        atomic_num = pt.GetAtomicNumber(symbol)
        if atomic_num <= 0:
            return None
        v = pt.GetDefaultValence(atomic_num)
        return v if v > 0 else None
    except Exception:
        return None


def _is_plausible_formula(elements: Dict[str, int]) -> bool:
    """Generic valence-based plausibility check for a candidate species.

    Sum of per-atom default valences must be at least ``2*(N-1)`` so the
    atoms can form a connected graph. Rejects H-overloaded species like
    ``*CH5`` / ``*NH4`` that blind element-arithmetic fallbacks can emit.
    Unknown-valence atoms (transition metals, exotic elements) bypass the
    check and are accepted — no chemistry-specific rules are hardcoded.
    """
    elements = {e: c for e, c in elements.items() if c and c > 0}
    if not elements:
        return False
    valences: Dict[str, int] = {}
    for e in elements:
        v = _element_default_valence(e)
        if v is None:
            return True
        valences[e] = v
    n_atoms = sum(elements.values())
    if n_atoms == 1:
        return True
    total_valence = sum(valences[e] * c for e, c in elements.items())
    return total_valence >= 2 * (n_atoms - 1)


def _filter_mcts_candidates(
    candidates: List[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """Drop candidates whose product formula fails the valence check.

    This is applied only inside the MCTS engine before candidates are added
    to the search tree or used for rollout sampling. Pruning-strategy
    enumeration is untouched.
    """
    kept: List[Dict[str, Any]] = []
    for c in candidates:
        prod_elements = c.get("product_elements")
        if not prod_elements:
            prod_elements = parse_species_elements(c.get("product_label", ""))
        if _is_plausible_formula(prod_elements):
            kept.append(c)
    return kept


@dataclass
class MCTSNode:
    species_label: str
    elements: Dict[str, int] = field(default_factory=dict)
    depth: int = 0
    parent: Optional["MCTSNode"] = None
    step_from_parent: Optional[Dict[str, Any]] = None

    children: List["MCTSNode"] = field(default_factory=list)
    untried_candidates: Optional[List[Dict[str, Any]]] = None

    visits: int = 0
    total_reward: float = 0.0
    is_terminal: bool = False  # reached target
    is_dead_end: bool = False  # no children possible

    @property
    def avg_reward(self) -> float:
        return self.total_reward / self.visits if self.visits else 0.0

    def ucb1(self, c: float) -> float:
        # Sign convention: rewards are -max(ΔG) in eV, so LOWER energy
        # paths have HIGHER avg_reward, and selection maximises UCB1.
        # With eV-scale rewards (ΔG_ads convention) the default c=1.4
        # exploration term is commensurate with the exploitation term.
        if self.visits == 0:
            return float("inf")
        parent_n = self.parent.visits if self.parent else 1
        exploit = self.avg_reward
        explore = c * math.sqrt(math.log(max(parent_n, 1)) / self.visits)
        return exploit + explore

    def is_fully_expanded(self) -> bool:
        return self.untried_candidates is not None and len(self.untried_candidates) == 0


class MCTSSearchEngine:
    """Monte Carlo Tree Search over catalytic state space.

    Each iteration:
      1. **Select** a leaf by descending UCB1 from the root.
      2. **Expand** one untried child if the leaf is non-terminal.
      3. **Simulate** a full rollout from that leaf until it hits the target
         (success) or a dead end (failure).
      4. **Backpropagate** the rollout reward up the tree. Successful
         rollouts also emit their full (root→rollout-tail) path as a
         candidate pathway.

    The ``rollout_depth`` parameter is a safety cap on rollout length, not a
    hard stop: without it, a cyclic or pathological graph could loop forever,
    but under normal conditions the rollout terminates on target/dead-end.
    """

    def __init__(
        self,
        generator: Optional[CandidateGeneratorRouter] = None,
        crn_graph: Optional[Dict[str, List[Dict[str, Any]]]] = None,
        energy_cache: Optional[EnergyCache] = None,
        energy_fn: Optional[Callable[[str], Optional[float]]] = None,
        step_dg_fn: Optional[Callable[[str, str, Dict[str, int], Dict[str, int]], Optional[float]]] = None,
        c_uct: float = 1.4,
        rollout_depth: int = 32,
        max_depth: int = 10,
        max_iterations: int = 100,
        max_state_evaluations: int = 40,
        rollout_bias_beta: float = 1.0,
        expansion_bias_beta: float = 0.5,
        distance_penalty_eV: float = 0.3,
        seed: Optional[int] = None,
    ):
        self.generator = generator or CandidateGeneratorRouter(enable_care=False)
        # Pre-built CRN adjacency dict. When set, ``_get_children`` does a
        # dict lookup instead of regenerating candidates via the router
        # every rollout step — the O(1) lookup is ~1000× faster than the
        # formula-merge loop over the full registry.
        self._crn_graph: Optional[Dict[str, List[Dict[str, Any]]]] = crn_graph
        self.energy_cache = energy_cache or EnergyCache()
        # (label) -> Optional[float]; must return the corrected adsorption
        # free energy ΔG_ads (shared eV-scale convention) or None when the
        # species is unevaluable. The energy_fn owns the convention for
        # gas-phase / composite labels too (e.g. ΔG_ads(X(g)) = 0 by
        # definition, composite = sum of components).
        self.energy_fn = energy_fn
        # When provided, reward is computed as -max(ΔG_step) along the path
        # (stoichiometry-corrected — preferred). Without it, reward falls
        # back to -max(ΔG_ads) over the path's evaluable nodes.
        self.step_dg_fn = step_dg_fn
        self.c_uct = c_uct
        self.rollout_depth = rollout_depth
        self.max_depth = max_depth
        self.max_iterations = max_iterations
        self.max_state_evaluations = max_state_evaluations
        # Target-aware sampling: at rollout/expansion time we pick candidates
        # weighted by ``exp(-beta * element_distance_to_target)`` so the
        # search preferentially moves toward the target composition. Set
        # beta=0 for pure random sampling. The same generic formula-space
        # distance works for any reaction / any surface.
        self.rollout_bias_beta = rollout_bias_beta
        self.expansion_bias_beta = expansion_bias_beta
        # Per-unit penalty (in eV-equivalent) for the element-count distance
        # between a rollout endpoint and the target composition. Without it,
        # rollouts that never reach the target all get identical zero reward
        # (when energies are unknown or flat) and UCB1 cannot distinguish
        # promising deep branches from unproductive side branches. Keep this
        # small relative to typical adsorption energies so the actual
        # energetic reward still dominates when available.
        self.distance_penalty_eV = distance_penalty_eV
        self._rng = random.Random(seed)

        self.root: Optional[MCTSNode] = None
        self._target_label_norm: str = ""
        self._target_elements: Dict[str, int] = {}
        self._eval_count = 0
        # Successful rollouts (those that reach the target) are remembered
        # here as full sequences of ``{species_label, elements, step_info}``
        # dicts so they can be emitted as candidate pathways even when no
        # tree node was expanded to the target.
        self._rollout_successes: List[List[Dict[str, Any]]] = []
        self._rollout_success_count: int = 0
        self._rollout_deadend_count: int = 0
        self._rollout_cap_count: int = 0

    def _get_children(
        self, species_label: str, elements: Dict[str, int],
    ) -> List[Dict[str, Any]]:
        """Return outgoing edges for a species.

        Uses the pre-built CRN adjacency dict when available (O(1)
        lookup), else falls back to calling the candidate generator
        on-demand. The on-demand fallback iterates over the full
        registry with valence checks and is ~1000× slower per call, so
        passing ``crn_graph`` to the engine is the right default.
        """
        if self._crn_graph is not None:
            adj = self._crn_graph.get(species_label, [])
            # Convert graph edge dicts to candidate dicts matching the
            # shape ``_filter_mcts_candidates`` and the rollout loop expect.
            return [
                {
                    "product_label": e["label"],
                    "product_elements": e["elements"],
                    **(e["step"] if isinstance(e.get("step"), dict) else {}),
                }
                for e in adj
            ]
        return self.generator.generate_candidates(
            species_label=species_label, elements=elements,
        )

    # ----- public -----

    def run(
        self,
        initial_label: str,
        target_label: str,
        initial_elements: Optional[Dict[str, int]] = None,
        target_elements: Optional[Dict[str, int]] = None,
    ) -> Dict[str, Any]:
        """Execute MCTS from ``initial_label`` toward ``target_label``.

        Returns ``{"paths": List[List[Node-dict]], "stats": {...}}``.
        """
        self.root = MCTSNode(
            species_label=initial_label,
            elements=initial_elements or parse_species_elements(initial_label),
            depth=0,
        )
        self._target_label_norm = self._strip_phase(target_label)
        self._target_elements = (
            target_elements or parse_species_elements(target_label)
        )
        self._eval_count = 0
        self._rollout_successes = []
        self._rollout_success_count = 0
        self._rollout_deadend_count = 0
        self._rollout_cap_count = 0

        terminal_nodes: List[MCTSNode] = []
        it = 0
        for it in range(self.max_iterations):
            if self._eval_count >= self.max_state_evaluations:
                logger.info("MCTS stop: eval budget %d reached at iter %d",
                            self.max_state_evaluations, it)
                break

            leaf = self._select(self.root)
            expanded = self._expand(leaf)
            node_to_simulate = expanded if expanded is not None else leaf
            reward = self._simulate(node_to_simulate)
            self._backpropagate(node_to_simulate, reward)

            if node_to_simulate.is_terminal and node_to_simulate not in terminal_nodes:
                terminal_nodes.append(node_to_simulate)

        paths = self._collect_paths(terminal_nodes, top_k=5)

        return {
            "paths": paths,
            "stats": {
                "iterations_run": it + 1 if self.max_iterations else 0,
                "state_evaluations": self._eval_count,
                "terminal_tree_nodes": len(terminal_nodes),
                "rollout_successes": self._rollout_success_count,
                "rollout_deadends": self._rollout_deadend_count,
                "rollout_cap_hits": self._rollout_cap_count,
                "cache": self.energy_cache.stats(),
            },
        }

    # ----- MCTS phases -----

    def _select(self, node: MCTSNode) -> MCTSNode:
        """Descend via UCB1 until a node that is not fully expanded or terminal."""
        cur = node
        while not cur.is_terminal and not cur.is_dead_end:
            if cur.untried_candidates is None:
                # Lazy candidate generation at first visit
                self._populate_candidates(cur)
            if cur.untried_candidates:
                return cur
            if not cur.children:
                cur.is_dead_end = True
                return cur
            cur = max(cur.children, key=lambda ch: ch.ucb1(self.c_uct))
        return cur

    def _expand(self, node: MCTSNode) -> Optional[MCTSNode]:
        if node.is_terminal or node.is_dead_end:
            return None
        if not node.untried_candidates:
            return None
        if node.depth >= self.max_depth:
            return None

        # Pick one candidate to expand; bias toward those whose product
        # formula is closer to the target composition.
        idx = self._pick_candidate_index(
            node.untried_candidates, self.expansion_bias_beta,
        )
        cand = node.untried_candidates.pop(idx)

        product_label = cand.get("product_label", "")
        product_elements = cand.get("product_elements") or parse_species_elements(product_label)

        child = MCTSNode(
            species_label=product_label,
            elements=product_elements,
            depth=node.depth + 1,
            parent=node,
            step_from_parent=cand,
        )
        child.is_terminal = self._matches_target(child)
        node.children.append(child)
        return child

    def _simulate(self, node: MCTSNode) -> float:
        """One rollout: generate a full root→target pathway and score it.

        Two sub-steps (both part of one MCTS simulation):

        **(a) Path generation.** Starting from ``node``, repeatedly pick the
        next intermediate by target-composition–biased sampling (softmax
        over element-count distance) until the target composition is hit.
        No energy evaluations are done here — the candidate generator is
        generative, so a path to any target composition always exists and
        the sampler is guaranteed to converge as long as some candidate
        reduces the distance. A safety cap (``rollout_depth``) guards
        against pathological cycles.

        **(b) Energy evaluation.** Once the full root→target path is in
        hand, evaluate every intermediate via the shared energy cache
        (cache hits are free; misses call ``energy_fn`` and charge the
        per-run state-evaluation budget). The reward is the negative of
        the worst (bottleneck) energy on the realised path — the MCTS
        "game score" for this simulation.

        Successful rollouts are stored in ``self._rollout_successes`` so
        they can be surfaced as candidate pathways even when no tree node
        was expanded all the way down to the target.
        """
        # --- (a) Generate the full path ---
        tree_chain = self._node_chain(node)
        visited: Set[str] = {
            self._strip_phase(n.species_label) for n in tree_chain
        }

        rollout_chain: List[Dict[str, Any]] = []
        cur_label = node.species_label
        cur_elements = node.elements

        reached_target = node.is_terminal
        steps = 0
        stuck = False
        while not reached_target and steps < self.rollout_depth:
            cands = self._get_children(cur_label, cur_elements)
            cands = _filter_mcts_candidates(cands)
            cands = [
                c for c in cands
                if self._strip_phase(c.get("product_label", "")) not in visited
            ]
            if not cands:
                stuck = True
                break

            idx = self._pick_candidate_index(cands, self.rollout_bias_beta)
            cand = cands[idx]
            cur_label = cand.get("product_label", "")
            cur_elements = (
                cand.get("product_elements")
                or parse_species_elements(cur_label)
            )
            visited.add(self._strip_phase(cur_label))
            rollout_chain.append({
                "species_label": cur_label,
                "elements": cur_elements,
                "step_info": cand,
            })
            steps += 1
            if self._matches_target_label(cur_label, cur_elements):
                reached_target = True

        if reached_target:
            self._rollout_success_count += 1
        elif stuck:
            self._rollout_deadend_count += 1
        else:
            self._rollout_cap_count += 1

        # --- (b) Evaluate energies on the assembled path ---
        # Each node on the path is an intermediate species; energy lookup is
        # cache-first so repeated intermediates across rollouts cost nothing.
        all_labels = [n.species_label for n in tree_chain] + [
            r["species_label"] for r in rollout_chain
        ]
        all_elements = [n.elements for n in tree_chain] + [
            r["elements"] for r in rollout_chain
        ]
        # Prefer ΔG_step bottleneck (stoichiometry-corrected, size-fair).
        # Fall back to the ΔG_ads bottleneck (same eV-scale convention)
        # when step_dg_fn is unavailable. Either way the reward is on the
        # eV scale — never an absolute UMA total.
        if self.step_dg_fn is not None:
            bottleneck = self._bottleneck_step_dG(all_labels, all_elements)
        else:
            bottleneck = self._bottleneck_dG_ads(all_labels)

        # Stitch and record the pathway for the rollout-successes listing.
        if reached_target:
            stitched: List[Dict[str, Any]] = []
            for n in tree_chain:
                stitched.append({
                    "species_label": n.species_label,
                    "elements": n.elements,
                    "step_info": n.step_from_parent or {},
                    "visits": n.visits,
                    "avg_reward": n.avg_reward,
                })
            for r in rollout_chain:
                stitched.append({
                    "species_label": r["species_label"],
                    "elements": r["elements"],
                    "step_info": r["step_info"],
                    "visits": 0,  # rollout-only, not a tree node
                    "avg_reward": 0.0,
                })
            self._rollout_successes.append(stitched)

        # Reward: negative of the worst step ΔG (or worst ΔG_ads) on the
        # full path — lower energy = higher reward. When NOTHING on the
        # path was evaluable (bottleneck = -inf) the path is unevaluable:
        # it carries no energetic signal and only the distance-to-target
        # penalty below differentiates it. For the rare rollouts that did
        # not reach target the penalty also steers UCB1 toward branches
        # whose trajectories end close to the goal.
        energy_reward = -bottleneck if bottleneck != float("-inf") else 0.0
        if not reached_target and self._target_elements and self.distance_penalty_eV > 0.0:
            end_distance = self._element_distance(cur_elements, self._target_elements)
            return energy_reward - self.distance_penalty_eV * end_distance
        return energy_reward

    def _backpropagate(self, node: MCTSNode, reward: float) -> None:
        cur: Optional[MCTSNode] = node
        while cur is not None:
            cur.visits += 1
            cur.total_reward += reward
            cur = cur.parent

    # ----- helpers -----

    def _populate_candidates(self, node: MCTSNode) -> None:
        if node.is_terminal or node.depth >= self.max_depth:
            node.untried_candidates = []
            return
        cands = self._get_children(node.species_label, node.elements)
        cands = _filter_mcts_candidates(cands)
        # Skip revisits along current root-to-node path
        path_norm = {self._strip_phase(l) for l in self._labels_from_root(node)}
        filtered = [
            c for c in cands
            if self._strip_phase(c.get("product_label", "")) not in path_norm
        ]
        node.untried_candidates = filtered

    def _lookup_energy(self, label: str) -> Optional[float]:
        """Cache-first ΔG_ads lookup. The convention for ALL labels —
        including gas-phase '(g)' and composite '+' ones — is owned by
        ``energy_fn`` (see ``__init__``); this method never special-cases
        them, so it cannot poison a shared cache with placeholder values.
        """
        hit, val = self.energy_cache.get(label)
        if hit:
            return val
        if self.energy_fn is None:
            self.energy_cache.put(label, None)
            return None
        val = self.energy_fn(label)
        self.energy_cache.put(label, val)
        self._eval_count += 1
        return val

    def _bottleneck_dG_ads(self, labels: List[str]) -> float:
        """Largest ΔG_ads along the path (eV scale, shared convention).

        Fallback reward basis when ``step_dg_fn`` is unavailable: the
        highest corrected adsorption free energy the path passes through.
        Unevaluable species (None) do not contribute. Returns -inf when no
        node on the path is evaluable (path unevaluable)."""
        worst = float("-inf")
        for lab in labels:
            e = self._lookup_energy(lab)
            if e is not None and e > worst:
                worst = e
        return worst

    def _bottleneck_step_dG(
        self,
        labels: List[str],
        elements_list: List[Dict[str, int]],
    ) -> float:
        """Largest ΔG_step (parent→child) along the path. Negative reward
        of this is the kinetically/thermodynamically meaningful signal:
        the highest energy uphill the path forces the system to climb."""
        worst = float("-inf")
        for i in range(len(labels) - 1):
            dG = self.step_dg_fn(
                labels[i], labels[i + 1],
                elements_list[i], elements_list[i + 1],
            )
            if dG is not None and dG > worst:
                worst = dG
        return worst

    def _labels_from_root(self, node: MCTSNode) -> List[str]:
        return [n.species_label for n in self._node_chain(node)]

    def _node_chain(self, node: MCTSNode) -> List[MCTSNode]:
        chain: List[MCTSNode] = []
        cur: Optional[MCTSNode] = node
        while cur is not None:
            chain.append(cur)
            cur = cur.parent
        chain.reverse()
        return chain

    def _matches_target(self, node: MCTSNode) -> bool:
        return self._matches_target_label(node.species_label, node.elements)

    def _matches_target_label(self, label: str, elements: Dict[str, int]) -> bool:
        if self._strip_phase(label) == self._target_label_norm:
            return True
        if self._target_elements and elements == self._target_elements:
            return True
        return False

    # ----- target-aware sampling helpers -----

    @staticmethod
    def _element_distance(
        current: Dict[str, int], target: Dict[str, int],
    ) -> int:
        """L1 distance between two element-count dicts.

        Generic (no chemistry specific to any reaction): the number of atoms
        that must be added or removed to transform ``current`` into
        ``target``. Minimum is 0 (identical composition).
        """
        if not target:
            return 0
        keys = set(current) | set(target)
        return sum(
            abs(current.get(e, 0) - target.get(e, 0)) for e in keys
        )

    def _pick_candidate_index(
        self,
        candidates: List[Dict[str, Any]],
        beta: float,
    ) -> int:
        """Pick an index from ``candidates`` weighted by distance-to-target.

        Weight for candidate i is ``exp(-beta * d_i)`` where ``d_i`` is the
        element-count distance from its product formula to the target. When
        ``beta <= 0`` or no target is set, falls back to uniform random —
        preserving pure MCTS behavior when target-awareness is disabled.
        """
        n = len(candidates)
        if n == 1:
            return 0
        if beta <= 0.0 or not self._target_elements:
            return self._rng.randrange(n)

        # Shift-safe softmax over distance to avoid overflow for deep targets.
        distances = []
        for c in candidates:
            pe = c.get("product_elements") or parse_species_elements(
                c.get("product_label", "")
            )
            distances.append(self._element_distance(pe, self._target_elements))
        d_min = min(distances)
        weights = [math.exp(-beta * (d - d_min)) for d in distances]
        total = sum(weights)
        if total <= 0.0:
            return self._rng.randrange(n)
        r = self._rng.random() * total
        acc = 0.0
        for i, w in enumerate(weights):
            acc += w
            if acc >= r:
                return i
        return n - 1

    @staticmethod
    def _strip_phase(label: str) -> str:
        s = label.strip().lower().replace(" ", "")
        return s.replace("*", "").replace("(g)", "").replace("(s)", "").replace("(l)", "")

    # ----- path extraction -----

    def _collect_paths(
        self, terminals: List[MCTSNode], top_k: int = 5,
    ) -> List[List[Dict[str, Any]]]:
        """Return top-K target-reaching paths.

        Includes paths found two ways:
          (a) Tree nodes that matched the target after expansion, and
          (b) Rollout simulations that reached the target (captured during
              ``_simulate`` into ``self._rollout_successes``).

        Deduplicated by label sequence. Ranked PRIMARILY by the energetic
        objective — highest reward, i.e. lowest pathway max-ΔG (rewards
        are -max(ΔG_step) / -max(ΔG_ads) in eV) — with path length only
        as a tiebreaker (shorter = simpler mechanism preferred).
        """
        candidates: List[Tuple[int, float, List[Dict[str, Any]]]] = []

        for t in terminals:
            chain = [
                {
                    "species_label": n.species_label,
                    "elements": n.elements,
                    "step_info": n.step_from_parent or {},
                    "visits": n.visits,
                    "avg_reward": n.avg_reward,
                }
                for n in self._node_chain(t)
            ]
            candidates.append((len(chain), t.avg_reward, chain))

        for chain in self._rollout_successes:
            # Score rollout-found paths by their tree prefix's avg_reward;
            # only the tree-portion has meaningful visit counts.
            tree_reward = 0.0
            for step in chain:
                if step.get("visits", 0) > 0:
                    tree_reward = step.get("avg_reward", 0.0)
            candidates.append((len(chain), tree_reward, chain))

        # Deduplicate by label sequence signature.
        seen: Set[Tuple[str, ...]] = set()
        unique: List[Tuple[int, float, List[Dict[str, Any]]]] = []
        for length, reward, chain in candidates:
            sig = tuple(s["species_label"] for s in chain)
            if sig in seen:
                continue
            seen.add(sig)
            unique.append((length, reward, chain))

        # Rank by energetic objective first (higher reward = lower pathway
        # max-ΔG); path length is only a tiebreaker.
        unique.sort(key=lambda x: (-x[1], x[0]))
        return [chain for _, _, chain in unique[:top_k]]


__all__ = ["MCTSNode", "MCTSSearchEngine"]
