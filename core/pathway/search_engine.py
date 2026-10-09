"""Mechanism search engine: best-first + beam search over catalytic states.

Expands candidate states layer by layer, evaluates free energies, prunes
siblings, and extracts complete pathways from root to goal.
"""

from __future__ import annotations

import logging
from collections import deque
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Set

from core.pathway.candidate_generators import (
    CandidateGeneratorRouter,
    parse_species_elements,
)
from core.pathway.free_energy_router import FreeEnergyRouter
from core.pathway.pruning_policy import PruningPolicy

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Search node
# ---------------------------------------------------------------------------

@dataclass
class PathSearchNode:
    """One node in the search tree."""
    state_id: str
    species_label: str
    elements: Dict[str, int] = field(default_factory=dict)
    depth: int = 0
    parent_id: Optional[str] = None
    step_from_parent: Optional[Dict[str, Any]] = None
    free_energy_eV: Optional[float] = None
    is_goal: bool = False
    is_pruned: bool = False
    prune_reason: str = ""


# ---------------------------------------------------------------------------
# Search engine
# ---------------------------------------------------------------------------

class MechanismSearchEngine:
    """Layer-by-layer beam search over catalytic state space."""

    def __init__(
        self,
        generator: Optional[CandidateGeneratorRouter] = None,
        energy_router: Optional[FreeEnergyRouter] = None,
        policy: Optional[PruningPolicy] = None,
    ):
        self.generator = generator or CandidateGeneratorRouter(enable_care=False)
        self.energy_router = energy_router or FreeEnergyRouter()
        self.policy = policy or PruningPolicy()

        # Search state
        self._nodes: Dict[str, PathSearchNode] = {}
        self._frontier: List[str] = []
        self._evaluations_used: int = 0
        self._goal_nodes: List[str] = []
        self._visited_labels: Set[str] = set()  # species_label dedup across depths

    # ----- public API -----

    def initialize_root(
        self,
        species_label: str,
        elements: Optional[Dict[str, int]] = None,
        free_energy_eV: Optional[float] = None,
    ) -> str:
        """Create the root node for the search."""
        if elements is None:
            elements = parse_species_elements(species_label)

        state_id = self._make_state_id(species_label, depth=0)
        root = PathSearchNode(
            state_id=state_id,
            species_label=species_label,
            elements=elements,
            depth=0,
            free_energy_eV=free_energy_eV,
        )
        self._nodes[state_id] = root
        self._frontier = [state_id]
        self._visited_labels.add(self._visited_key(species_label))
        return state_id

    def expand_node(self, node_id: str) -> List[str]:
        """Expand a single node into its children."""
        node = self._nodes.get(node_id)
        if node is None:
            return []

        candidates = self.generator.generate_candidates(
            species_label=node.species_label,
            elements=node.elements,
        )

        child_ids: List[str] = []
        for cand in candidates:
            product_label = cand.get("product_label", "")
            product_elements = cand.get("product_elements", {})
            if not product_elements:
                product_elements = parse_species_elements(product_label)

            child_id = self._make_state_id(product_label, depth=node.depth + 1)
            # Skip duplicates (by state_id AND by normalized label across
            # depths). The visited key PRESERVES phase markers: '*CH4'
            # (adsorbed) and 'CH4(g)' (desorbed) are distinct states, so
            # visiting one must not block creating the other.
            #
            # Known limitation (first-parent-wins): a species reached via
            # one parent is never re-created for a different parent, even
            # if the alternative partial path is lower in free energy —
            # energies are only evaluated after expansion, so a safe
            # energy-based revisit rule is not available at this point.
            # Alternative routes to the same intermediate are therefore
            # represented only by the first one discovered.
            normalized = self._visited_key(product_label)
            if child_id in self._nodes or normalized in self._visited_labels:
                continue
            self._visited_labels.add(normalized)

            child = PathSearchNode(
                state_id=child_id,
                species_label=product_label,
                elements=product_elements,
                depth=node.depth + 1,
                parent_id=node_id,
                step_from_parent=cand,
            )
            self._nodes[child_id] = child
            child_ids.append(child_id)

        return child_ids

    def evaluate_frontier(
        self,
        surface_path: Optional[str] = None,
        backend: Optional[str] = None,
        fixed_site: Optional[list] = None,
    ) -> int:
        """Evaluate corrected adsorption free energies (ΔG_ads) for all
        frontier nodes.

        ``free_energy_eV`` is either ΔG_ads on the shared convention (see
        ``core.pathway.free_energy_router`` module docstring) or None when
        the state cannot be evaluated. None means INCOMPARABLE — pruning
        keeps such nodes and never ranks them against numeric values. No
        0.0 placeholders are ever stored.
        """
        count = 0
        for nid in self._frontier:
            node = self._nodes.get(nid)
            if node is None or node.free_energy_eV is not None:
                continue

            # Gas-phase and composite labels cannot be evaluated by the
            # single-adsorbate UMA path here: leave free_energy_eV=None
            # (incomparable; kept by pruning), never a fabricated value.
            label = node.species_label
            if "(g)" in label or "+" in label:
                logger.debug(
                    "evaluate_frontier: %s is gas/composite — "
                    "ΔG_ads unavailable at this layer, kept as incomparable.",
                    label,
                )
                continue

            extra_kwargs: Dict[str, Any] = {}
            if fixed_site is not None:
                extra_kwargs["fixed_site"] = fixed_site

            est = self.energy_router.estimate_state_free_energy(
                label,
                backend=backend,
                surface_path=surface_path,
                **extra_kwargs,
            )
            node.free_energy_eV = est.get("free_energy_eV")
            self._evaluations_used += 1
            count += 1

        # Release GPU memory after each frontier evaluation batch
        try:
            import torch
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
        except Exception:
            pass

        return count

    def prune_frontier(self) -> Dict[str, Any]:
        """Apply sibling comparison + beam width to the current frontier."""
        # Group siblings by parent
        parent_groups: Dict[Optional[str], List[str]] = {}
        for nid in self._frontier:
            node = self._nodes[nid]
            parent = node.parent_id
            parent_groups.setdefault(parent, []).append(nid)

        all_kept: List[str] = []
        all_pruned: List[str] = []
        decisions: Dict[str, str] = {}

        for parent, children in parent_groups.items():
            estimates = []
            for cid in children:
                node = self._nodes[cid]
                estimates.append({
                    "state_id": cid,
                    "free_energy_eV": node.free_energy_eV,
                })

            result = self.policy.compare_siblings(estimates)
            kept = result["kept"] + result["marginal"]
            pruned = result["pruned"]

            for pid in pruned:
                pnode = self._nodes.get(pid)
                if pnode:
                    pnode.is_pruned = True
                    pnode.prune_reason = result["decisions"].get(pid, "pruned")

            all_kept.extend(kept)
            all_pruned.extend(pruned)
            decisions.update(result["decisions"])

        # Apply beam width
        kept_nodes = [{"state_id": k, "free_energy_eV": self._nodes[k].free_energy_eV} for k in all_kept]
        selected, dropped = self.policy.select_frontier(kept_nodes)
        selected_ids = [s["state_id"] for s in selected]

        for d in dropped:
            did = d["state_id"]
            dnode = self._nodes.get(did)
            if dnode:
                dnode.is_pruned = True
                dnode.prune_reason = "dropped by beam width"
            all_pruned.append(did)

        self._frontier = selected_ids
        return {
            "kept": selected_ids,
            "pruned": all_pruned,
            "decisions": decisions,
        }

    def is_goal_state(
        self,
        node_id: str,
        target_label: str,
        target_elements: Optional[Dict[str, int]] = None,
    ) -> bool:
        """Check if a node matches the target state.

        Handles phase-agnostic matching: ``*CH4`` matches ``CH4(g)`` because
        they share the same element composition.
        """
        node = self._nodes.get(node_id)
        if node is None:
            return False

        # Label match (normalized — strips *, (g), (s), whitespace)
        if self._normalize_label(node.species_label) == self._normalize_label(target_label):
            return True

        # Element match (phase-agnostic)
        if target_elements and node.elements == target_elements:
            return True

        # Phase-stripped label match: *CH4 vs CH4(g) → both normalize to "ch4"
        if self._strip_phase(node.species_label) == self._strip_phase(target_label):
            return True

        return False

    def backtrack_path(self, goal_id: str) -> List[PathSearchNode]:
        """Trace back from a goal node to the root, returning the path."""
        path: List[PathSearchNode] = []
        current_id: Optional[str] = goal_id
        while current_id is not None:
            node = self._nodes.get(current_id)
            if node is None:
                break
            path.append(node)
            current_id = node.parent_id
        path.reverse()
        return path

    def search_until_stop(
        self,
        target_label: str,
        target_elements: Optional[Dict[str, int]] = None,
        surface_path: Optional[str] = None,
        backend: Optional[str] = None,
        fixed_site: Optional[list] = None,
    ) -> Dict[str, Any]:
        """Run the full search loop.

        Returns a dict with ``pathways`` (list of node sequences),
        ``stats``, and ``all_nodes``.
        """
        if target_elements is None:
            target_elements = parse_species_elements(target_label)

        depth = 0
        while True:
            # Check stopping condition
            should_stop, reason = self.policy.should_stop(depth, self._evaluations_used)
            if should_stop:
                logger.info("Search stopped: %s", reason)
                break

            if not self._frontier:
                logger.info("Search stopped: empty frontier")
                break

            # Expand all frontier nodes
            new_children: List[str] = []
            for nid in list(self._frontier):
                children = self.expand_node(nid)
                new_children.extend(children)

            if not new_children:
                logger.info("Search stopped: no new candidates at depth %d", depth)
                break

            # Check for goal states among children
            for cid in new_children:
                if self.is_goal_state(cid, target_label, target_elements):
                    self._nodes[cid].is_goal = True
                    self._goal_nodes.append(cid)

            # Move frontier to children
            self._frontier = new_children

            # Evaluate and prune
            self.evaluate_frontier(surface_path=surface_path, backend=backend, fixed_site=fixed_site)
            self.prune_frontier()

            depth += 1

        # Extract pathways from goal nodes
        pathways = []
        for gid in self._goal_nodes:
            path = self.backtrack_path(gid)
            pathways.append(path)

        # If no goal reached, extract best frontier paths as partial candidates
        if not pathways and self._frontier:
            frontier_nodes = [
                (nid, self._nodes[nid])
                for nid in self._frontier
                if nid in self._nodes
            ]
            # Sort by energy (lowest first)
            frontier_nodes.sort(
                key=lambda x: x[1].free_energy_eV if x[1].free_energy_eV is not None else float("inf"),
            )
            # Take top-3 best frontier paths
            for nid, _ in frontier_nodes[:3]:
                path = self.backtrack_path(nid)
                if len(path) > 1:
                    pathways.append(path)
            if pathways:
                logger.info(
                    "No goal reached; returning %d best partial paths (depth %d)",
                    len(pathways), depth,
                )

        return {
            "pathways": pathways,
            "stats": {
                "depth_reached": depth,
                "evaluations_used": self._evaluations_used,
                "total_nodes": len(self._nodes),
                "goal_nodes_found": len(self._goal_nodes),
                "frontier_size": len(self._frontier),
            },
            "all_nodes": self._nodes,
        }

    # ----- graph enumeration for systematic search -----

    def build_crn_graph(
        self,
        max_depth: int = 10,
        max_graph_nodes: int = 2000,
        target_label: Optional[str] = None,
        target_elements: Optional[Dict[str, int]] = None,
    ) -> Dict[str, List[Dict[str, Any]]]:
        """Stage-1 CRN construction only: return an adjacency dict
        ``{species_label: [{label, elements, step}, ...]}``.

        Fast BFS graph build with no energy evaluation and no path
        extraction. Meant for Stage-2 algorithms (pruning or MCTS) that
        navigate the pre-built graph via O(1) neighbour lookup instead
        of calling the generator on-demand.

        When ``target_label`` (and optionally ``target_elements``) is
        given, the graph is post-filtered to only species that lie on
        at least one path from root to target — species that are dead
        ends or irrelevant to the target are removed. The filter is
        element-agnostic: it just checks label / element equivalence
        and graph reachability.
        """
        graph: Dict[str, List[Dict[str, Any]]] = {}
        if not self._frontier:
            return graph
        root_node = self._nodes[self._frontier[0]]
        expand_queue: deque = deque([
            (root_node.species_label, root_node.elements, 0),
        ])
        # Phase-preserving visited key: '*X' and 'X(g)' are distinct states.
        visited: Set[str] = {self._visited_key(root_node.species_label)}

        while expand_queue:
            label, elements, depth = expand_queue.popleft()
            if depth >= max_depth:
                continue
            children = self.generator.generate_candidates(
                species_label=label, elements=elements,
            )
            adj = []
            for c in children:
                prod = c.get("product_label", "")
                prod_elem = c.get("product_elements", {})
                if not prod_elem:
                    prod_elem = parse_species_elements(prod)
                adj.append({"label": prod, "elements": prod_elem, "step": c})

                norm = self._visited_key(prod)
                if norm not in visited:
                    if len(visited) >= max_graph_nodes:
                        continue
                    visited.add(norm)
                    expand_queue.append((prod, prod_elem, depth + 1))
            graph[label] = adj

        logger.info(
            "CRN graph built: %d species, %d outgoing edge lists",
            len(visited), len(graph),
        )

        if target_label is not None:
            graph = self._filter_graph_to_target(
                graph, root_node.species_label,
                target_label, target_elements,
            )
            logger.info(
                "CRN graph filtered to target: %d species, %d outgoing edge lists",
                len(set(graph.keys()) | {c["label"] for adj in graph.values() for c in adj}),
                len(graph),
            )
        return graph

    def _filter_graph_to_target(
        self,
        graph: Dict[str, List[Dict[str, Any]]],
        root_label: str,
        target_label: str,
        target_elements: Optional[Dict[str, int]] = None,
    ) -> Dict[str, List[Dict[str, Any]]]:
        """Return subgraph containing only species that appear on some
        path from ``root_label`` to any target-matching species.

        Algorithm:
          1. Forward reachable from root (given: all species in graph).
          2. Backward reachable from target: reverse BFS via inverse
             edges collected from the adjacency dict.
          3. Keep intersection; drop species / edges outside it.
        """
        target_norm = self._strip_phase(target_label)

        def _is_target(lbl: str, elem: Dict[str, int]) -> bool:
            if self._strip_phase(lbl) == target_norm:
                return True
            if target_elements and elem == target_elements:
                return True
            return False

        # Collect every species' element dict in a single pass (used for
        # target matching below). This replaces an earlier O(N²·E)
        # double loop over species × graph-edges.
        species_elements: Dict[str, Dict[str, int]] = {}
        for label, adj in graph.items():
            for e in adj:
                species_elements.setdefault(e["label"], e["elements"])
        species_elements.setdefault(root_label, {})
        all_species: Set[str] = set(species_elements.keys())
        all_species.update(graph.keys())

        # Find target species: match by normalized label, or by element dict.
        target_species: Set[str] = set()
        for sp in all_species:
            if self._strip_phase(sp) == target_norm:
                target_species.add(sp)
            elif target_elements and species_elements.get(sp) == target_elements:
                target_species.add(sp)
        if not target_species:
            logger.info("CRN target not reachable in graph; returning full graph.")
            return graph

        # Build reverse adjacency: for each (u, v) edge, reverse[v].add(u)
        reverse: Dict[str, Set[str]] = {}
        for u, adj in graph.items():
            for e in adj:
                reverse.setdefault(e["label"], set()).add(u)

        # Reverse BFS from target(s)
        backward_reachable: Set[str] = set(target_species)
        queue: deque = deque(target_species)
        while queue:
            node = queue.popleft()
            for pred in reverse.get(node, ()):
                if pred not in backward_reachable:
                    backward_reachable.add(pred)
                    queue.append(pred)

        # Forward reachable: BFS from root (only consider species in graph)
        forward_reachable: Set[str] = {root_label}
        queue = deque([root_label])
        while queue:
            node = queue.popleft()
            for e in graph.get(node, ()):
                if e["label"] not in forward_reachable:
                    forward_reachable.add(e["label"])
                    queue.append(e["label"])

        keep = forward_reachable & backward_reachable
        filtered: Dict[str, List[Dict[str, Any]]] = {}
        for u, adj in graph.items():
            if u not in keep:
                continue
            kept_adj = [e for e in adj if e["label"] in keep]
            if kept_adj:
                filtered[u] = kept_adj
        return filtered

    def compute_dist_to_target(
        self,
        graph: Dict[str, List[Dict[str, Any]]],
        target_label: str,
        target_elements: Optional[Dict[str, int]] = None,
    ) -> Dict[str, int]:
        """Reverse-BFS distance (in CRN hops) from every node to target.

        Returns a dict ``{label: hops_to_target}``. Nodes that cannot
        reach target are absent. Used by pruning's distance-bucket beam
        and min-dist sibling exemption to guarantee target-progress
        nodes survive energy pruning.
        """
        # Build reverse adjacency + collect each node's element dict.
        reverse_adj: Dict[str, Set[str]] = {}
        node_elements: Dict[str, Dict[str, int]] = {}
        for u, edges in graph.items():
            reverse_adj.setdefault(u, set())
            for e in edges:
                reverse_adj.setdefault(e["label"], set()).add(u)
                node_elements.setdefault(e["label"], e["elements"])
        all_labels: Set[str] = set(reverse_adj.keys()) | set(graph.keys())

        target_norm = self._strip_phase(target_label)
        target_labels: Set[str] = set()
        for lbl in all_labels:
            if self._strip_phase(lbl) == target_norm:
                target_labels.add(lbl)
            elif target_elements and node_elements.get(lbl) == target_elements:
                target_labels.add(lbl)

        dist: Dict[str, int] = {t: 0 for t in target_labels}
        queue: deque = deque(target_labels)
        while queue:
            node = queue.popleft()
            d = dist[node]
            for pred in reverse_adj.get(node, ()):
                if pred not in dist:
                    dist[pred] = d + 1
                    queue.append(pred)
        return dist

    def enumerate_all_paths_bfs(
        self,
        target_label: str,
        target_elements: Optional[Dict[str, int]] = None,
        max_depth: int = 10,
        max_paths: int = 50,
        max_graph_nodes: int = 2000,
        max_dfs_steps: int = 200000,
    ) -> List[List[Dict[str, Any]]]:
        """CRN construction + bounded path extraction (Stage 1 of UniMech).

        Pure graph-level work with NO energy evaluation. Thin wrapper:
        delegates graph build to ``build_crn_graph`` and path extraction
        to ``extract_paths_from_graph`` so both implementations share one
        source of truth.
        """
        if not self._frontier:
            return []
        root_node = self._nodes[self._frontier[0]]
        if target_elements is None:
            target_elements = parse_species_elements(target_label)

        graph = self.build_crn_graph(
            max_depth=max_depth,
            max_graph_nodes=max_graph_nodes,
            target_label=None,  # no target filter — Stage 1 wants the full graph
            target_elements=None,
        )
        return self.extract_paths_from_graph(
            graph,
            root_label=root_node.species_label,
            root_elements=root_node.elements,
            target_label=target_label,
            target_elements=target_elements,
            max_depth=max_depth,
            max_paths=max_paths,
            max_dfs_steps=max_dfs_steps,
        )

    def extract_paths_from_graph(
        self,
        graph: Dict[str, List[Dict[str, Any]]],
        root_label: str,
        root_elements: Dict[str, int],
        target_label: str,
        target_elements: Optional[Dict[str, int]] = None,
        max_depth: int = 10,
        max_paths: int = 50,
        max_dfs_steps: int = 200000,
    ) -> List[List[Dict[str, Any]]]:
        """Extract initial→target paths from a pre-built CRN adjacency dict.

        Pure DFS path extraction — no graph construction here, no energy
        evaluation. Used by Stage 1 to derive the path list for the
        ``quick`` strategy from the same CRN dict that ``pruning`` /
        ``mcts`` traverse, so all three see identical search input.
        """
        if target_elements is None:
            target_elements = parse_species_elements(target_label)

        all_paths: List[List[Dict[str, Any]]] = []
        dfs_steps = [0]
        dfs_budget_hit = [False]

        def _elem_distance(cur: Dict[str, int]) -> int:
            if not target_elements:
                return 0
            keys = set(cur) | set(target_elements)
            return sum(
                abs(cur.get(e, 0) - target_elements.get(e, 0)) for e in keys
            )

        def _dfs(current: str, current_elem: Dict[str, int],
                 path: List[Dict[str, Any]], visited: Set[str]):
            if len(all_paths) >= max_paths or dfs_budget_hit[0]:
                return
            if len(path) > max_depth + 1:
                return
            dfs_steps[0] += 1
            if dfs_steps[0] > max_dfs_steps:
                dfs_budget_hit[0] = True
                return

            if self._strip_phase(current) == self._strip_phase(target_label):
                all_paths.append(list(path))
                return
            if target_elements and current_elem == target_elements:
                all_paths.append(list(path))
                return

            remaining = max_depth - (len(path) - 1)
            for child in graph.get(current, []):
                if dfs_budget_hit[0] or len(all_paths) >= max_paths:
                    return
                child_label = child["label"]
                # Phase-preserving key: '*X' → 'X(g)' (desorption) is a
                # legitimate step within one path.
                child_norm = self._visited_key(child_label)
                if child_norm in visited:
                    continue
                child_elem = child["elements"]
                child_dist = _elem_distance(child_elem)
                if child_dist > max(remaining - 1, 0) * 6:
                    continue
                visited.add(child_norm)
                path.append({
                    "species_label": child_label,
                    "elements": child_elem,
                    "step_info": child.get("step", {}),
                })
                _dfs(child_label, child_elem, path, visited)
                path.pop()
                visited.remove(child_norm)

        root_entry = {
            "species_label": root_label,
            "elements": root_elements,
            "step_info": {},
        }
        _dfs(
            root_label, root_elements,
            [root_entry],
            {self._visited_key(root_label)},
        )

        if dfs_budget_hit[0]:
            logger.info(
                "extract_paths_from_graph: DFS budget %d exhausted; "
                "returning %d paths.",
                max_dfs_steps, len(all_paths),
            )
        logger.info(
            "extract_paths_from_graph: %d DFS steps, %d paths to %s",
            dfs_steps[0], len(all_paths), target_label,
        )
        return all_paths

    # ----- internal helpers -----

    @staticmethod
    def _make_state_id(label: str, depth: int) -> str:
        """Create a unique state ID."""
        safe = label.replace("*", "ads_").replace("(", "_").replace(")", "_").replace("+", "_")
        return f"d{depth}_{safe}"

    @staticmethod
    def _normalize_label(label: str) -> str:
        return label.strip().lower().replace(" ", "")

    @staticmethod
    def _strip_phase(label: str) -> str:
        """Remove phase markers: '*CO' → 'co', 'CH4(g)' → 'ch4'.

        Only for phase-AGNOSTIC matching (e.g. goal tests). Do NOT use as
        a visited/dedup key — see ``_visited_key``.
        """
        s = label.strip().lower().replace(" ", "")
        s = s.replace("*", "").replace("(g)", "").replace("(s)", "").replace("(l)", "")
        return s

    @staticmethod
    def _visited_key(label: str) -> str:
        """Dedup key that PRESERVES phase markers: '*CH4' and 'CH4(g)' are
        distinct thermodynamic states (adsorbed vs desorbed) and must both
        be creatable during the search."""
        return label.strip().lower().replace(" ", "")


__all__ = ["PathSearchNode", "MechanismSearchEngine"]
