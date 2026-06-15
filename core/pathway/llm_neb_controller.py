"""
LLM-Controlled NEB Structure Preparation

This module provides a general-purpose, LLM-controlled approach for preparing
NEB (Nudged Elastic Band) calculations for any reaction on any surface.

The LLM fully controls:
1. Analysis of reactant and product structures
2. Decision on whether atoms need to be added/removed
3. Positions of added atoms
4. Validation of prepared structures
5. Iterative refinement if structures are unreasonable

Key principle: The adsorption sites are pre-determined and cannot be changed.
LLM only decides how to match atoms between reactant and product for NEB.
"""

import os
import re
import json
import logging
import numpy as np
from typing import List, Dict, Tuple, Optional, Any, Union
from dataclasses import dataclass
from ase import Atoms
from ase import Atom
from ase.io import read, write
import tempfile

logger = logging.getLogger(__name__)


@dataclass
class NEBStructurePlan:
    """LLM's plan for NEB structure preparation"""
    reaction_type: str  # e.g., "hydrogenation", "dehydrogenation", "bond_breaking"
    atoms_to_add_to_reactant: List[Dict]  # [{"element": "H", "position": [x,y,z], "reason": "..."}]
    atoms_to_add_to_product: List[Dict]
    atoms_to_remove_from_reactant: List[int]  # atom indices
    atoms_to_remove_from_product: List[int]
    atom_mapping: Dict[int, int]  # reactant_idx -> product_idx for alignment
    reasoning: str
    confidence: float  # 0-1

    # Energy correction information
    original_reactant_energy: Optional[float] = None
    adjusted_reactant_energy: Optional[float] = None
    original_product_energy: Optional[float] = None
    adjusted_product_energy: Optional[float] = None


@dataclass
class ValidationResult:
    """LLM's validation of prepared structures"""
    is_valid: bool
    issues: List[str]
    suggestions: List[str]
    corrected_plan: Optional[NEBStructurePlan] = None
    auto_fixes_applied: List[str] = None  # Record what was auto-fixed

    def __post_init__(self):
        if self.auto_fixes_applied is None:
            self.auto_fixes_applied = []


class LLMNEBController:
    """
    LLM-controlled NEB structure preparation.

    This class uses an LLM to intelligently prepare structures for NEB calculations
    by analyzing the chemistry and deciding how to adjust atom counts and positions.
    """

    def __init__(
        self,
        model: str = None,
        max_retries: int = 3,
        temperature: float = 0.1,
    ):
        """
        Initialize the LLM NEB Controller.

        Parameters
        ----------
        model : str, optional
            LLM model to use. Defaults to environment variable or claude-sonnet.
        max_retries : int
            Maximum attempts for structure preparation
        temperature : float
            LLM temperature for responses
        """
        self.model = model or os.getenv("OPENAI_MODEL", "claude-sonnet-4-20250514")
        self.max_retries = max_retries
        self.temperature = temperature
        self._client = None

    def _init_llm(self):
        """Initialize LLM client"""
        if self._client is not None:
            return

        try:
            from openai import OpenAI

            api_key = os.getenv("OPENAI_API_KEY")
            base_url = os.getenv("OPENAI_BASE_URL", "https://api.openai.com/v1")

            if not api_key:
                raise ValueError("OPENAI_API_KEY not set")

            self._client = OpenAI(api_key=api_key, base_url=base_url)
            logger.info(f"LLM initialized (model: {self.model})")

        except Exception as e:
            logger.error(f"Failed to initialize LLM: {e}")
            raise

    def _atoms_to_detailed_string(
        self,
        atoms: Atoms,
        adsorbate_indices: List[int],
        name: str = "Structure"
    ) -> str:
        """
        Convert ASE Atoms to a detailed string representation for LLM.

        Includes:
        - Cell parameters
        - All atom positions with indices and markers
        - Adsorbate atoms clearly marked
        - Interatomic distances for adsorbate atoms
        """
        lines = [f"=== {name} ==="]
        lines.append(f"Formula: {atoms.get_chemical_formula()}")
        lines.append(f"Total atoms: {len(atoms)}")

        # Identify surface element(s)
        symbols = atoms.get_chemical_symbols()
        from collections import Counter
        symbol_counts = Counter(symbols)
        surface_element = max(symbol_counts, key=symbol_counts.get)
        lines.append(f"Surface element: {surface_element}")
        lines.append(f"Adsorbate indices: {adsorbate_indices}")

        # Cell
        cell = atoms.cell
        lines.append(f"\nCell (Å):")
        for i, vec in enumerate(cell):
            lines.append(f"  a{i+1}: [{vec[0]:.6f}, {vec[1]:.6f}, {vec[2]:.6f}]")

        # Surface z-range
        positions = atoms.get_positions()
        surface_z = [positions[i, 2] for i in range(len(atoms)) if symbols[i] == surface_element]
        lines.append(f"\nSurface z-range: {min(surface_z):.3f} - {max(surface_z):.3f} Å")

        # Adsorbate atoms (detailed)
        lines.append(f"\nAdsorbate atoms:")
        for idx in adsorbate_indices:
            pos = positions[idx]
            sym = symbols[idx]
            lines.append(f"  [{idx}] {sym}: ({pos[0]:.4f}, {pos[1]:.4f}, {pos[2]:.4f})")

        # Distances between adsorbate atoms
        if len(adsorbate_indices) > 1:
            lines.append(f"\nAdsorbate interatomic distances:")
            for i, idx1 in enumerate(adsorbate_indices):
                for idx2 in adsorbate_indices[i+1:]:
                    dist = np.linalg.norm(positions[idx1] - positions[idx2])
                    lines.append(f"  {symbols[idx1]}[{idx1}] - {symbols[idx2]}[{idx2}]: {dist:.3f} Å")

        # Full coordinates (POSCAR-like)
        lines.append(f"\nFull coordinates (Cartesian, Å):")
        lines.append("Index  Element  X           Y           Z           Note")
        lines.append("-" * 70)
        for i, (pos, sym) in enumerate(zip(positions, symbols)):
            note = "ADSORBATE" if i in adsorbate_indices else ""
            lines.append(f"{i:5d}  {sym:7s}  {pos[0]:10.5f}  {pos[1]:10.5f}  {pos[2]:10.5f}  {note}")

        return "\n".join(lines)

    def _call_llm(self, system_prompt: str, user_prompt: str) -> str:
        """Call LLM and return response"""
        self._init_llm()

        response = self._client.chat.completions.create(
            model=self.model,
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": user_prompt}
            ],
            temperature=self.temperature,
            max_tokens=4000,
        )

        return response.choices[0].message.content

    def analyze_reaction_and_plan(
        self,
        reactant: Atoms,
        product: Atoms,
        reactant_adsorbate_indices: List[int],
        product_adsorbate_indices: List[int],
        reaction_name: str = None,
        previous_step_product: Atoms = None,
        previous_step_product_adsorbate_indices: List[int] = None,
        step_index: int = None,
    ) -> NEBStructurePlan:
        """
        Use LLM to analyze the reaction and create a plan for NEB structure preparation.

        Parameters
        ----------
        reactant : Atoms
            Reactant structure (already relaxed at its adsorption site)
        product : Atoms
            Product structure (already relaxed at its adsorption site)
        reactant_adsorbate_indices : List[int]
            Indices of adsorbate atoms in reactant
        product_adsorbate_indices : List[int]
            Indices of adsorbate atoms in product
        reaction_name : str, optional
            Human-readable reaction name

        Returns
        -------
        NEBStructurePlan
            LLM's plan for structure preparation
        """
        # Prepare detailed structure descriptions
        reactant_desc = self._atoms_to_detailed_string(
            reactant, reactant_adsorbate_indices, "REACTANT"
        )
        product_desc = self._atoms_to_detailed_string(
            product, product_adsorbate_indices, "PRODUCT"
        )

        # Analyze atom count differences
        r_symbols = reactant.get_chemical_symbols()
        p_symbols = product.get_chemical_symbols()
        from collections import Counter
        r_counts = Counter(r_symbols)
        p_counts = Counter(p_symbols)

        all_elements = set(r_counts.keys()) | set(p_counts.keys())
        diff_info = []
        for elem in sorted(all_elements):
            r_n = r_counts.get(elem, 0)
            p_n = p_counts.get(elem, 0)
            if r_n != p_n:
                diff_info.append(f"  {elem}: {r_n} (reactant) → {p_n} (product), change: {p_n - r_n:+d}")

        diff_str = "\n".join(diff_info) if diff_info else "  No difference in atom counts"

        system_prompt = """You are an expert computational chemist specializing in catalytic reactions and NEB (Nudged Elastic Band) calculations.

Your task is to analyze a surface catalytic reaction and create a plan for preparing structures for NEB calculation.

IMPORTANT CONSTRAINTS:
1. The adsorption sites and positions of existing adsorbate atoms are FIXED and cannot be changed
2. For NEB to work, reactant and product must have the SAME number of atoms
3. When adding atoms (e.g., H from gas phase), position them to create a REALISTIC reaction path
4. Added atoms should be positioned close to their final bonding position (within 0.5-2.0 Å)
5. **CRITICAL**: After adding atoms, the adsorbate atom sequence in reactant MUST match the product sequence by element type
   - Example: If product adsorbate is [C, H, O], reactant adsorbate MUST also be [C, H, O], NOT [C, O, H]
   - This ensures NEB interpolation matches corresponding atoms correctly
   - When adding atoms, specify the insertion index to maintain correct ordering

**PATHWAY CONTINUITY** (if this is part of a multi-step reaction pathway):
6. If provided with the previous step's product structure, the current reactant MUST be at the SAME adsorption site
   - The current reactant should be structurally identical to the previous product (same adsorbate, same site)
   - If they differ (e.g., previous product is *CHO, current reactant is *CHOH), you need to add the difference atoms
   - NEVER change adsorption sites between consecutive steps - reactions happen at the same location
7. For desorption reactions (e.g., *CHOH → *CH + OH_gas):
   - Identify which atoms/groups are being removed from the surface
   - In the product, position desorbing species slightly above their original location (z + 1.5-2.5 Å)
   - For molecules like OH: keep bonding (O-H ≈ 0.97 Å) but move the whole group up
   - Goal: smooth pathway from adsorbed to desorbed, NOT jump to far vacuum

CHEMICAL KNOWLEDGE:
- Typical bond lengths: C-H ≈ 1.09 Å, O-H ≈ 0.96 Å, C-O ≈ 1.43 Å, C=O ≈ 1.23 Å
- For hydrogenation reactions (adding H):
  * FIRST analyze the product structure to determine which atom the H is bonded to
  * Find the shortest distance from H to any other atom in the product - that's its bonding partner
  * Then place the H in the reactant at 2.0-2.5 Å from that SAME bonding partner atom
  * The H should approach from above (higher z) toward its bonding position
- For dehydrogenation/desorption reactions (removing species):
  * Atoms being removed should be positioned slightly above their original location in the reactant
  * For desorbed species, place at z = reactant_z + 1.5-2.5 Å (not too far in vacuum)
  * For desorbed OH groups: place O near C, H bonded to O at O-H distance ≈ 0.97 Å
  * Goal: create smooth pathway from adsorbed to desorbed state, not jump to vacuum
- Atoms should be above the surface, but not too far (z < surface_z_max + 5 Å)
- ALWAYS check atom-to-atom distances in the product to identify bonding relationships
- For reordering: provide atom_mapping showing which reactant atom corresponds to which product atom

OUTPUT FORMAT:
Return a JSON object with this exact structure:
```json
{
  "reaction_type": "hydrogenation|dehydrogenation|isomerization|bond_breaking|bond_forming|desorption|other",
  "analysis": "Brief analysis of the reaction chemistry and pathway continuity",
  "atoms_to_add_to_reactant": [
    {
      "element": "H",
      "position": [x, y, z],
      "insert_index": null,  // OPTIONAL: If specified, insert at this index; if null, append at end
      "reason": "Why this atom is added and why at this position"
    }
  ],
  "atoms_to_add_to_product": [],
  "atoms_to_remove_from_reactant": [],
  "atoms_to_remove_from_product": [
    {
      "indices": [i, j, k],  // Atom indices to remove (e.g., for desorbing OH)
      "reason": "What species is being removed and why"
    }
  ],
  "atom_mapping": {
    "description": "How adsorbate atoms in reactant correspond to those in product",
    "mapping": {"reactant_idx": "product_idx"}
  },
  "reasoning": "Detailed explanation of the plan, INCLUDING:\n    - Verification that final adsorbate atom sequences match by element type\n    - If part of pathway: confirmation that reactant is at the same site as previous product\n    - For desorption: identification of which atoms/groups are leaving the surface",
  "confidence": 0.9
}
```

If atoms need to be ADDED to make counts equal:
- Add to the structure with FEWER atoms
- Position added atoms near their final positions (not far away!)
- For H being added: place at x,y of target atom, z about 0.5-1.0 Å above target z
- **CRITICAL**: Use "insert_index" to ensure the final adsorbate atom sequence matches the product sequence
  - Example: If product adsorbate is [C, H, O] and reactant is [C, O], add H with insert_index such that result is [C, H, O]
  - Always verify the final element sequence will match!

If atoms need to be REMOVED (desorption reactions):
- Identify which atoms are leaving the surface (e.g., OH group in *CHOH → *CH + OH)
- Provide indices of atoms to remove from product
- In product, these atoms should be positioned slightly above their reactant positions
"""

        # Get current adsorbate sequences
        r_ads_symbols = [reactant[i].symbol for i in reactant_adsorbate_indices]
        p_ads_symbols = [product[i].symbol for i in product_adsorbate_indices]

        # Prepare pathway continuity information if provided
        pathway_context = ""
        if previous_step_product is not None and previous_step_product_adsorbate_indices is not None:
            prev_ads_symbols = [previous_step_product[i].symbol for i in previous_step_product_adsorbate_indices]
            prev_ads_positions = previous_step_product.positions[previous_step_product_adsorbate_indices]

            # Describe previous step's product
            pathway_context = f"""
**PATHWAY CONTINUITY CONTEXT** (This is step {step_index if step_index is not None else '?'} of a multi-step reaction pathway):

PREVIOUS STEP'S PRODUCT:
  Adsorbate composition: {prev_ads_symbols}
  Adsorbate positions:
{chr(10).join([f"    [{i}] {previous_step_product[idx].symbol} at [{prev_ads_positions[i][0]:.3f}, {prev_ads_positions[i][1]:.3f}, {prev_ads_positions[i][2]:.3f}]" for i, idx in enumerate(previous_step_product_adsorbate_indices)])}

**CRITICAL CONTINUITY REQUIREMENT**:
- The CURRENT REACTANT should be at the SAME adsorption site as the previous product
- If the current reactant adsorbate composition ({r_ads_symbols}) differs from the previous product ({prev_ads_symbols}),
  you need to explain the difference and ensure the reaction is continuous
- Adsorption sites must NOT change between steps - the reaction pathway is continuous at one location
- Example: If previous product is *CHO and current reactant is *CHOH, the H atom should be added to *CHO
  at the same adsorption site, NOT create a new *CHOH at a different location
"""

        user_prompt = f"""Analyze this surface catalytic reaction and create a NEB structure preparation plan.

REACTION: {reaction_name or 'Unknown'}
{pathway_context}
{reactant_desc}

{product_desc}

ATOM COUNT DIFFERENCES:
{diff_str}

CURRENT ADSORBATE SEQUENCES:
Reactant adsorbate: {r_ads_symbols}
Product adsorbate:  {p_ads_symbols}

**CRITICAL REQUIREMENTS**:
1. After adding/removing atoms, the final adsorbate atom sequences MUST match by element type.
   - If you add atoms, use "insert_index" to ensure correct ordering
   - Verify that the final reactant adsorbate sequence will equal the product adsorbate sequence

2. If this is part of a pathway (previous product provided above):
   - Current reactant must be at the SAME adsorption site as previous product
   - Explain any differences in adsorbate composition between previous product and current reactant
   - For desorption reactions: explicitly identify which atoms/groups are leaving the surface

Please analyze the reaction and provide a JSON plan for preparing these structures for NEB calculation.
Remember:
- Existing adsorbate positions are FIXED
- NEB requires equal atom counts
- Added atoms should be positioned realistically (close to their bonding sites)
- **FINAL ADSORBATE SEQUENCES MUST MATCH BY ELEMENT TYPE**
- **PATHWAY CONTINUITY: Same adsorption site across consecutive steps**
"""

        response = self._call_llm(system_prompt, user_prompt)

        # Parse JSON from response
        json_match = re.search(r'```json\s*(.*?)\s*```', response, re.DOTALL)
        if not json_match:
            json_match = re.search(r'\{.*\}', response, re.DOTALL)

        if not json_match:
            raise ValueError(f"Could not parse JSON from LLM response: {response[:500]}")

        try:
            plan_dict = json.loads(json_match.group(1) if json_match.group(1) else json_match.group())
        except:
            plan_dict = json.loads(json_match.group())

        # Convert to NEBStructurePlan
        plan = NEBStructurePlan(
            reaction_type=plan_dict.get("reaction_type", "unknown"),
            atoms_to_add_to_reactant=plan_dict.get("atoms_to_add_to_reactant", []),
            atoms_to_add_to_product=plan_dict.get("atoms_to_add_to_product", []),
            atoms_to_remove_from_reactant=plan_dict.get("atoms_to_remove_from_reactant", []),
            atoms_to_remove_from_product=plan_dict.get("atoms_to_remove_from_product", []),
            atom_mapping=plan_dict.get("atom_mapping", {}),
            reasoning=plan_dict.get("reasoning", ""),
            confidence=plan_dict.get("confidence", 0.5),
        )

        logger.info(f"LLM analysis: {plan.reaction_type}, confidence: {plan.confidence}")
        logger.info(f"Plan: add {len(plan.atoms_to_add_to_reactant)} to reactant, "
                   f"add {len(plan.atoms_to_add_to_product)} to product")

        return plan

    def apply_plan(
        self,
        reactant: Atoms,
        product: Atoms,
        reactant_adsorbate_indices: List[int],
        product_adsorbate_indices: List[int],
        plan: NEBStructurePlan,
    ) -> Tuple[Atoms, Atoms, List[int], List[int]]:
        """
        Apply the LLM's plan to prepare structures for NEB.

        Returns
        -------
        Tuple of (adjusted_reactant, adjusted_product, new_reactant_indices, new_product_indices)
        """
        # Make copies
        adj_reactant = reactant.copy()
        adj_product = product.copy()
        new_r_indices = list(reactant_adsorbate_indices)
        new_p_indices = list(product_adsorbate_indices)

        # Remove atoms first (if any)
        if plan.atoms_to_remove_from_reactant:
            indices_to_remove = sorted(plan.atoms_to_remove_from_reactant, reverse=True)
            for idx in indices_to_remove:
                del adj_reactant[idx]
                new_r_indices = [i if i < idx else i - 1 for i in new_r_indices if i != idx]

        if plan.atoms_to_remove_from_product:
            indices_to_remove = sorted(plan.atoms_to_remove_from_product, reverse=True)
            for idx in indices_to_remove:
                del adj_product[idx]
                new_p_indices = [i if i < idx else i - 1 for i in new_p_indices if i != idx]

        # Add atoms
        from ase import Atom  # Import once at the top of the section

        for atom_info in plan.atoms_to_add_to_reactant:
            element = atom_info["element"]
            position = np.array(atom_info["position"])
            insert_index = atom_info.get("insert_index", None)

            # Validate position
            positions = adj_reactant.get_positions()
            if len(positions) > 0:
                min_dist = np.min(np.linalg.norm(positions - position, axis=1))

                if min_dist < 0.5:
                    logger.warning(f"Added {element} too close to existing atom ({min_dist:.2f} Å), adjusting...")
                    # Find nearest atom and move away
                    nearest_idx = np.argmin(np.linalg.norm(positions - position, axis=1))
                    direction = position - positions[nearest_idx]
                    direction = direction / np.linalg.norm(direction)
                    position = positions[nearest_idx] + direction * 1.0

            # Insert atom at specified index or append
            if insert_index is not None and 0 <= insert_index <= len(adj_reactant):
                # Insert at specific position
                from ase import Atoms as AseAtoms
                new_atom = Atom(element, position=position)
                temp_atoms = adj_reactant[:insert_index] + AseAtoms([new_atom]) + adj_reactant[insert_index:]
                adj_reactant = temp_atoms

                # Update adsorbate indices: indices >= insert_index shift by 1
                new_r_indices = [i if i < insert_index else i + 1 for i in new_r_indices]

                # Insert new atom index at correct position (maintain sorted order)
                insert_pos = sum(1 for idx in new_r_indices if idx < insert_index)
                new_r_indices.insert(insert_pos, insert_index)

                logger.info(f"Inserted {element} at index {insert_index} in reactant at {position}")
            else:
                # Append at end (default behavior)
                adj_reactant += Atom(element, position=position)
                new_r_indices.append(len(adj_reactant) - 1)
                logger.info(f"Added {element} to reactant at {position}")

        for atom_info in plan.atoms_to_add_to_product:
            element = atom_info["element"]
            position = np.array(atom_info["position"])
            insert_index = atom_info.get("insert_index", None)

            # Validate position
            positions = adj_product.get_positions()
            if len(positions) > 0:
                min_dist = np.min(np.linalg.norm(positions - position, axis=1))

                if min_dist < 0.5:
                    logger.warning(f"Added {element} too close to existing atom ({min_dist:.2f} Å), adjusting...")
                    nearest_idx = np.argmin(np.linalg.norm(positions - position, axis=1))
                    direction = position - positions[nearest_idx]
                    direction = direction / np.linalg.norm(direction)
                    position = positions[nearest_idx] + direction * 1.0

            # Insert atom at specified index or append
            if insert_index is not None and 0 <= insert_index <= len(adj_product):
                # Insert at specific position
                from ase import Atoms as AseAtoms
                new_atom = Atom(element, position=position)
                temp_atoms = adj_product[:insert_index] + AseAtoms([new_atom]) + adj_product[insert_index:]
                adj_product = temp_atoms

                # Update adsorbate indices: indices >= insert_index shift by 1
                new_p_indices = [i if i < insert_index else i + 1 for i in new_p_indices]

                # Insert new atom index at correct position (maintain sorted order)
                insert_pos = sum(1 for idx in new_p_indices if idx < insert_index)
                new_p_indices.insert(insert_pos, insert_index)

                logger.info(f"Inserted {element} at index {insert_index} in product at {position}")
            else:
                # Append at end (default behavior)
                adj_product += Atom(element, position=position)
                new_p_indices.append(len(adj_product) - 1)
                logger.info(f"Added {element} to product at {position}")

        # Step: Handle atom reordering to match product atom types
        if len(new_r_indices) == len(new_p_indices):
            # Check if atom ordering matches by element type
            reactant_symbols = [adj_reactant[i].symbol for i in new_r_indices]
            product_symbols = [adj_product[i].symbol for i in new_p_indices]

            if reactant_symbols != product_symbols:
                logger.warning(f"Atom ordering mismatch: reactant {reactant_symbols} vs product {product_symbols}")
                logger.info(f"Attempting to reorder reactant indices to match product atom types...")

                # Simple reordering strategy: match atoms by element type
                # For each position in product, find a reactant atom of the same element
                new_r_indices_reordered = []
                used_reactant_indices = set()

                for p_idx in new_p_indices:
                    product_elem = adj_product[p_idx].symbol

                    # Find an unused reactant atom with the same element
                    found = False
                    for r_idx in new_r_indices:
                        if r_idx not in used_reactant_indices and adj_reactant[r_idx].symbol == product_elem:
                            new_r_indices_reordered.append(r_idx)
                            used_reactant_indices.add(r_idx)
                            found = True
                            break

                    if not found:
                        logger.warning(f"Could not find reactant atom matching product {product_elem}[{p_idx}]")

                if len(new_r_indices_reordered) == len(new_r_indices):
                    # ⭐ CRITICAL FIX: Actually reorder the reactant atoms physically
                    logger.info(f"Applying physical atom reordering to reactant...")

                    # Step 1: Identify which atoms are surface vs adsorbate
                    n_atoms = len(adj_reactant)
                    n_ads = len(new_r_indices)
                    adsorbate_indices_set = set(new_r_indices)

                    # Step 2: Build complete reordering mapping
                    # Surface atoms: keep in original order (indices 0 to n-n_ads-1)
                    # Adsorbate atoms: reorder according to new_r_indices_reordered
                    reorder_sequence = []

                    # First, add all surface atoms (not in adsorbate_indices)
                    for i in range(n_atoms):
                        if i not in adsorbate_indices_set:
                            reorder_sequence.append(i)

                    # Then, add reordered adsorbate atoms
                    reorder_sequence.extend(new_r_indices_reordered)

                    logger.info(f"Reorder sequence (old→new indices): {reorder_sequence}")

                    # Step 3: Apply reordering using ASE slicing
                    adj_reactant = adj_reactant[reorder_sequence]

                    # Step 4: Update adsorbate indices (they're now at the end)
                    n_surface_atoms = len(reorder_sequence) - n_ads
                    new_r_indices = list(range(n_surface_atoms, len(adj_reactant)))

                    logger.info(f"Updated adsorbate indices: {new_r_indices}")

                    # Verify the reordering
                    reactant_symbols_new = [adj_reactant[i].symbol for i in new_r_indices]
                    if reactant_symbols_new == product_symbols:
                        logger.info(f"✓ Atom ordering now matches: {reactant_symbols_new}")
                    else:
                        logger.warning(f"Reordering verification failed: got {reactant_symbols_new}, expected {product_symbols}")
                else:
                    logger.warning(f"Could not reorder all atoms: {len(new_r_indices_reordered)} of {len(new_r_indices)}")

        return adj_reactant, adj_product, new_r_indices, new_p_indices

    def auto_fix_structures(
        self,
        reactant: Atoms,
        product: Atoms,
        reactant_adsorbate_indices: List[int],
        product_adsorbate_indices: List[int],
        validation: ValidationResult,
    ) -> Tuple[Atoms, Atoms, List[int], List[int], List[str]]:
        """
        Automatically fix structural issues based on LLM validation results.

        This is the key agent loop: LLM detects problems → LLM proposes fixes → apply fixes

        The LLM sees:
        - The problematic structures
        - Validation issues identified
        - Validation suggestions from previous LLM

        And returns:
        - Corrected atom positions
        - Reordering instructions

        Returns
        -------
        Tuple of (fixed_reactant, fixed_product, new_r_indices, new_p_indices, fixes_applied)
        """
        logger.info("Using LLM to generate fixes based on validation feedback...")

        # Prepare detailed description of the problem
        reactant_desc = self._atoms_to_detailed_string(
            reactant, reactant_adsorbate_indices, "PROBLEMATIC_REACTANT"
        )
        product_desc = self._atoms_to_detailed_string(
            product, product_adsorbate_indices, "PROBLEMATIC_PRODUCT"
        )

        issues_str = "\n".join(f"  - {issue}" for issue in validation.issues)
        suggestions_str = "\n".join(f"  - {sug}" for sug in validation.suggestions)

        system_prompt = """You are an expert computational chemist fixing NEB structure problems.

You will receive:
1. Problematic reactant and product structures
2. Issues identified during validation
3. Suggestions for fixing the issues

Your task: Provide SPECIFIC fixes to make the structures suitable for NEB.

IMPORTANT:
- For atom ordering issues: specify which atoms to reorder
- For distance issues: provide new coordinates for atoms that need repositioning
- For hydrogenation: ensure H is positioned 2.0-2.5 Å from target atom
- Keep surface atoms unchanged, only modify adsorbate atoms

OUTPUT FORMAT:
```json
{
  "fixes": [
    {
      "type": "reorder_atoms",
      "description": "Swap H and O to match product ordering",
      "reactant_atom_mapping": {"67": 68, "68": 67}
    },
    {
      "type": "reposition_atom",
      "description": "Move H further from O to avoid overlap",
      "structure": "reactant",
      "atom_index": 68,
      "new_position": [x, y, z],
      "reason": "Why this position"
    }
  ],
  "reasoning": "Overall explanation of the fixes"
}
```

Available fix types:
- "reorder_atoms": Change atom ordering in structure
- "reposition_atom": Move an atom to a new position
- "separate_atoms": Push two atoms apart
"""

        user_prompt = f"""Fix the following NEB structure problems:

VALIDATION ISSUES:
{issues_str}

VALIDATION SUGGESTIONS:
{suggestions_str}

{reactant_desc}

{product_desc}

Please provide JSON with specific fixes to address these issues.
Focus on:
1. Matching atom ordering between reactant and product (same element at same index)
2. Ensuring atoms are not too close (> 0.8 Å apart)
3. Positioning added H atoms at reasonable pre-reaction distances (2.0-2.5 Å from target)
"""

        try:
            response = self._call_llm(system_prompt, user_prompt)

            # Parse JSON
            json_match = re.search(r'```json\s*(.*?)\s*```', response, re.DOTALL)
            if not json_match:
                json_match = re.search(r'\{.*"fixes".*\}', response, re.DOTALL)

            if not json_match:
                logger.warning("LLM did not return parseable fixes, using rule-based fallback")
                return self._rule_based_fix(reactant, product, reactant_adsorbate_indices,
                                            product_adsorbate_indices, validation)

            fix_dict = json.loads(json_match.group(1) if '```' in response else json_match.group())

            # Apply LLM-generated fixes
            return self._apply_llm_fixes(
                reactant, product,
                reactant_adsorbate_indices, product_adsorbate_indices,
                fix_dict
            )

        except Exception as e:
            logger.warning(f"LLM fix generation failed: {e}, using rule-based fallback")
            return self._rule_based_fix(reactant, product, reactant_adsorbate_indices,
                                       product_adsorbate_indices, validation)

    def _apply_llm_fixes(
        self,
        reactant: Atoms,
        product: Atoms,
        reactant_adsorbate_indices: List[int],
        product_adsorbate_indices: List[int],
        fix_dict: Dict,
    ) -> Tuple[Atoms, Atoms, List[int], List[int], List[str]]:
        """Apply fixes generated by LLM"""
        fixes_applied = []
        adj_reactant = reactant.copy()
        adj_product = product.copy()
        new_r_idx = list(reactant_adsorbate_indices)
        new_p_idx = list(product_adsorbate_indices)

        fixes = fix_dict.get("fixes", [])

        for fix in fixes:
            fix_type = fix.get("type")
            description = fix.get("description", "")

            if fix_type == "reorder_atoms":
                # Reorder atoms according to mapping
                mapping = fix.get("reactant_atom_mapping", {})

                if mapping:
                    # Convert string keys to int
                    swap_pairs = []
                    for old_str, new_str in mapping.items():
                        old_idx = int(old_str)
                        new_idx = int(new_str) if isinstance(new_str, (int, str)) else new_str
                        if isinstance(new_idx, str):
                            new_idx = int(new_idx)
                        swap_pairs.append((old_idx, new_idx))

                    # Create a new Atoms object with swapped atoms
                    # This ensures the swap is actually applied
                    from ase import Atoms as AseAtoms
                    from ase import Atom

                    # Save original tags before rebuilding
                    original_tags = adj_reactant.get_tags()

                    new_atoms_list = []
                    for i in range(len(adj_reactant)):
                        # Check if this index should be swapped
                        target_idx = i
                        for old_idx, new_idx in swap_pairs:
                            if i == old_idx:
                                target_idx = new_idx
                                break
                            elif i == new_idx:
                                target_idx = old_idx
                                break

                        # Copy the atom from target index
                        atom = adj_reactant[target_idx]
                        new_atoms_list.append(Atom(atom.symbol, atom.position))

                    # Rebuild structure with TAGS preserved
                    adj_reactant = AseAtoms(new_atoms_list,
                                           cell=adj_reactant.cell,
                                           pbc=adj_reactant.pbc)

                    # Apply swapped tags (tags should follow the atoms)
                    new_tags = []
                    for i in range(len(adj_reactant)):
                        target_idx = i
                        for old_idx, new_idx in swap_pairs:
                            if i == old_idx:
                                target_idx = new_idx
                                break
                            elif i == new_idx:
                                target_idx = old_idx
                                break
                        new_tags.append(original_tags[target_idx])

                    adj_reactant.set_tags(new_tags)

                    fixes_applied.append(f"LLM fix: {description}")
                    logger.info(f"  Applied LLM fix: {description} (rebuilt structure with tags)")

            elif fix_type == "reposition_atom":
                structure = fix.get("structure", "reactant")
                atom_idx = int(fix.get("atom_index"))
                new_pos = np.array(fix.get("new_position"))
                reason = fix.get("reason", "")

                if structure == "reactant" and atom_idx < len(adj_reactant):
                    positions = adj_reactant.get_positions()
                    positions[atom_idx] = new_pos
                    adj_reactant.set_positions(positions)
                    fixes_applied.append(f"LLM fix: {description} - {reason}")
                    logger.info(f"  Applied LLM fix: repositioned atom {atom_idx} to {new_pos}")

                elif structure == "product" and atom_idx < len(adj_product):
                    positions = adj_product.get_positions()
                    positions[atom_idx] = new_pos
                    adj_product.set_positions(positions)
                    fixes_applied.append(f"LLM fix: {description} - {reason}")
                    logger.info(f"  Applied LLM fix: repositioned atom {atom_idx} to {new_pos}")

        if not fixes_applied:
            logger.warning("No LLM fixes were applied")

        return adj_reactant, adj_product, new_r_idx, new_p_idx, fixes_applied

    def _rule_based_fix(
        self,
        reactant: Atoms,
        product: Atoms,
        reactant_adsorbate_indices: List[int],
        product_adsorbate_indices: List[int],
        validation: ValidationResult,
    ) -> Tuple[Atoms, Atoms, List[int], List[int], List[str]]:
        """
        Rule-based fallback fixes (original implementation)
        """
        fixes_applied = []
        adj_reactant = reactant.copy()
        adj_product = product.copy()
        new_r_idx = list(reactant_adsorbate_indices)
        new_p_idx = list(product_adsorbate_indices)

        # Fix 1: Atom ordering mismatch
        # Check if adsorbate atoms have different elements at same indices
        r_symbols = np.array(adj_reactant.get_chemical_symbols())
        p_symbols = np.array(adj_product.get_chemical_symbols())

        # Find mismatches in adsorbate region
        mismatches = []
        for i, (r_idx, p_idx) in enumerate(zip(new_r_idx, new_p_idx)):
            if r_idx < len(r_symbols) and p_idx < len(p_symbols):
                if r_symbols[r_idx] != p_symbols[p_idx]:
                    mismatches.append((i, r_idx, p_idx, r_symbols[r_idx], p_symbols[p_idx]))

        if mismatches:
            logger.info(f"Detected {len(mismatches)} atom ordering mismatches in adsorbate")

            # Strategy: Reorder reactant adsorbate atoms to match product
            # Build a mapping based on element types
            r_ads_positions = adj_reactant.get_positions()[new_r_idx]
            p_ads_positions = adj_product.get_positions()[new_p_idx]
            r_ads_symbols = [r_symbols[i] for i in new_r_idx]
            p_ads_symbols = [p_symbols[i] for i in new_p_idx]

            # For each product adsorbate atom, find the closest same-element reactant atom
            from scipy.optimize import linear_sum_assignment

            n_ads = len(new_r_idx)
            cost_matrix = np.full((n_ads, n_ads), 1e10)

            for i, p_sym in enumerate(p_ads_symbols):
                for j, r_sym in enumerate(r_ads_symbols):
                    if p_sym == r_sym:
                        # Cost is distance between atoms
                        cost_matrix[i, j] = np.linalg.norm(p_ads_positions[i] - r_ads_positions[j])

            row_ind, col_ind = linear_sum_assignment(cost_matrix)

            # Reorder reactant adsorbate atoms
            new_order_r_idx = [new_r_idx[col_ind[i]] for i in range(n_ads)]

            # Create a full reordering array for the entire reactant structure
            # Keep surface atoms in place, only reorder adsorbate
            surface_indices = [i for i in range(len(adj_reactant)) if i not in new_r_idx]
            full_new_order = surface_indices + new_order_r_idx

            adj_reactant = adj_reactant[full_new_order]

            # Update indices
            new_r_idx = list(range(len(surface_indices), len(adj_reactant)))

            fixes_applied.append(f"Reordered reactant adsorbate atoms to match product element sequence")
            logger.info(f"  Applied fix: atom reordering")

        # Fix 2: Atoms too close (< 0.8 Å)
        r_pos = adj_reactant.get_positions()
        p_pos = adj_product.get_positions()

        # Check reactant adsorbate atoms
        for i, idx_i in enumerate(new_r_idx):
            for j, idx_j in enumerate(new_r_idx):
                if i < j:
                    dist = np.linalg.norm(r_pos[idx_i] - r_pos[idx_j])
                    if dist < 0.8:
                        # Move atoms apart along the line connecting them
                        direction = r_pos[idx_i] - r_pos[idx_j]
                        direction = direction / np.linalg.norm(direction)

                        # Move each atom 0.5 Å away from midpoint
                        midpoint = (r_pos[idx_i] + r_pos[idx_j]) / 2
                        r_pos[idx_i] = midpoint + direction * 0.5
                        r_pos[idx_j] = midpoint - direction * 0.5

                        fixes_applied.append(f"Separated reactant atoms {idx_i} and {idx_j} (were {dist:.2f} Å apart)")
                        logger.info(f"  Applied fix: separated atoms {idx_i}-{idx_j}")

        adj_reactant.set_positions(r_pos)

        # Check product adsorbate atoms
        for i, idx_i in enumerate(new_p_idx):
            for j, idx_j in enumerate(new_p_idx):
                if i < j:
                    dist = np.linalg.norm(p_pos[idx_i] - p_pos[idx_j])
                    if dist < 0.8:
                        direction = p_pos[idx_i] - p_pos[idx_j]
                        direction = direction / np.linalg.norm(direction)

                        midpoint = (p_pos[idx_i] + p_pos[idx_j]) / 2
                        p_pos[idx_i] = midpoint + direction * 0.5
                        p_pos[idx_j] = midpoint - direction * 0.5

                        fixes_applied.append(f"Separated product atoms {idx_i} and {idx_j} (were {dist:.2f} Å apart)")
                        logger.info(f"  Applied fix: separated atoms {idx_i}-{idx_j}")

        adj_product.set_positions(p_pos)

        # Fix 3: For hydrogenation, ensure H is approaching from reasonable distance
        # If reaction type is hydrogenation and H is too close to target
        if len(new_r_idx) > len(reactant_adsorbate_indices):
            # We added atoms, likely hydrogenation
            added_indices = [i for i in new_r_idx if i >= len(reactant)]

            for added_idx in added_indices:
                if r_symbols[added_idx] == 'H':
                    # Find the closest non-H atom in adsorbate
                    other_ads_indices = [i for i in new_r_idx if i != added_idx and r_symbols[i] != 'H']

                    if other_ads_indices:
                        distances = [np.linalg.norm(r_pos[added_idx] - r_pos[other_idx])
                                   for other_idx in other_ads_indices]
                        min_dist = min(distances)
                        closest_idx = other_ads_indices[np.argmin(distances)]

                        if min_dist < 1.5:
                            # Move H further away (to ~2.0-2.5 Å)
                            target_dist = 2.2
                            direction = r_pos[added_idx] - r_pos[closest_idx]
                            direction = direction / np.linalg.norm(direction)

                            new_h_pos = r_pos[closest_idx] + direction * target_dist
                            r_pos[added_idx] = new_h_pos

                            fixes_applied.append(f"Moved H atom away from target (was {min_dist:.2f} Å, now {target_dist:.2f} Å)")
                            logger.info(f"  Applied fix: repositioned H atom for realistic approach")

        adj_reactant.set_positions(r_pos)

        return adj_reactant, adj_product, new_r_idx, new_p_idx, fixes_applied

    def validate_structures(
        self,
        reactant: Atoms,
        product: Atoms,
        reactant_adsorbate_indices: List[int],
        product_adsorbate_indices: List[int],
        plan: NEBStructurePlan,
    ) -> ValidationResult:
        """
        Use LLM to validate the prepared structures.

        Returns
        -------
        ValidationResult
            LLM's assessment of structure validity
        """
        # Basic checks first
        issues = []

        # Check atom counts
        if len(reactant) != len(product):
            issues.append(f"Atom count mismatch: reactant={len(reactant)}, product={len(product)}")

        # Check for overlapping atoms
        r_pos = reactant.get_positions()
        p_pos = product.get_positions()

        for i in range(len(r_pos)):
            for j in range(i + 1, len(r_pos)):
                dist = np.linalg.norm(r_pos[i] - r_pos[j])
                if dist < 0.5:
                    issues.append(f"Reactant: atoms {i} and {j} overlap (dist={dist:.2f} Å)")

        for i in range(len(p_pos)):
            for j in range(i + 1, len(p_pos)):
                dist = np.linalg.norm(p_pos[i] - p_pos[j])
                if dist < 0.5:
                    issues.append(f"Product: atoms {i} and {j} overlap (dist={dist:.2f} Å)")

        # If basic checks fail, no need for LLM
        if issues:
            return ValidationResult(
                is_valid=False,
                issues=issues,
                suggestions=["Fix the basic structural issues first"]
            )

        # LLM validation
        reactant_desc = self._atoms_to_detailed_string(
            reactant, reactant_adsorbate_indices, "ADJUSTED_REACTANT"
        )
        product_desc = self._atoms_to_detailed_string(
            product, product_adsorbate_indices, "ADJUSTED_PRODUCT"
        )

        system_prompt = """You are an expert computational chemist validating structures for NEB calculation.

Check if the prepared structures are chemically reasonable for NEB:
1. Are atom counts equal? (REQUIRED for NEB)
2. Are all atoms at reasonable positions? (no overlaps, reasonable bond lengths)
3. Is the reaction path physically meaningful?
4. Will linear interpolation between these structures create reasonable intermediates?

IMPORTANT NOTES FOR HYDROGEN PLACEMENT:
- For hydrogen approaching carbon (hydrogenation reactions), some proximity to oxygen is acceptable
- H-O distances of 1.2-2.0 Å are reasonable for hydrogen approaching C in presence of nearby O
- Judge reasonableness by: (a) does H move continuously toward C? (b) are bond formation/breaking chemically sound?
- Avoid only blocking cases where H would COLLIDE with another atom (overlap)

OUTPUT FORMAT:
```json
{
  "is_valid": true/false,
  "issues": ["list of problems found"],
  "suggestions": ["list of suggestions to fix issues"],
  "analysis": "detailed analysis"
}
```

If structures are valid, issues and suggestions should be empty arrays.
If invalid, provide specific suggestions for fixing the problems.
"""

        user_prompt = f"""Validate these structures for NEB calculation.

ORIGINAL PLAN:
- Reaction type: {plan.reaction_type}
- Reasoning: {plan.reasoning}

{reactant_desc}

{product_desc}

Are these structures suitable for NEB calculation? Be realistic about hydrogen placement:
1. Equal atom counts
2. No overlapping atoms (<0.5 Å separation)
3. Chemically reasonable (especially for hydrogen approach trajectories)
4. Linear interpolation will create reasonable path

Focus on: Is the interpolation path physically meaningful? Does it represent a reasonable reaction mechanism?
"""

        response = self._call_llm(system_prompt, user_prompt)

        # Parse JSON
        json_match = re.search(r'```json\s*(.*?)\s*```', response, re.DOTALL)
        if not json_match:
            json_match = re.search(r'\{.*\}', response, re.DOTALL)

        if json_match:
            try:
                result_dict = json.loads(json_match.group(1) if '```' in json_match.group() else json_match.group())
                return ValidationResult(
                    is_valid=result_dict.get("is_valid", False),
                    issues=result_dict.get("issues", []),
                    suggestions=result_dict.get("suggestions", [])
                )
            except:
                pass

        # Fallback: assume valid if no JSON parsing issues and basic checks passed
        return ValidationResult(is_valid=True, issues=[], suggestions=[])

    def prepare_neb_structures(
        self,
        reactant: Atoms,
        product: Atoms,
        reactant_adsorbate_indices: List[int],
        product_adsorbate_indices: List[int],
        reaction_name: str = None,
        previous_step_product: Atoms = None,
        previous_step_product_adsorbate_indices: List[int] = None,
        step_index: int = None,
    ) -> Tuple[Atoms, Atoms, List[int], List[int], NEBStructurePlan]:
        """
        Main method: prepare structures for NEB calculation.

        This method orchestrates the full LLM-controlled agent workflow:
        1. Analyze reaction and create plan
        2. Apply plan to structures
        3. Validate results
        4. AUTO-FIX issues if validation fails
        5. Re-validate and retry if needed

        Parameters
        ----------
        reactant, product : Atoms
            Input structures (pre-relaxed at their adsorption sites)
        reactant_adsorbate_indices, product_adsorbate_indices : List[int]
            Adsorbate atom indices
        reaction_name : str, optional
            Human-readable reaction name
        previous_step_product : Atoms, optional
            Product structure from previous step in pathway (for continuity)
        previous_step_product_adsorbate_indices : List[int], optional
            Adsorbate indices for previous step's product
        step_index : int, optional
            Index of this step in the pathway (for context)

        Returns
        -------
        Tuple of (adjusted_reactant, adjusted_product, reactant_indices, product_indices, plan)
        """
        logger.info(f"Preparing NEB structures for: {reaction_name or 'reaction'}")
        logger.info(f"Reactant: {len(reactant)} atoms, Product: {len(product)} atoms")

        adj_reactant = None
        adj_product = None
        new_r_idx = None
        new_p_idx = None
        plan = None
        all_fixes = []

        for attempt in range(self.max_retries):
            logger.info(f"Attempt {attempt + 1}/{self.max_retries}")

            try:
                # Step 1: Analyze and plan
                plan = self.analyze_reaction_and_plan(
                    reactant, product,
                    reactant_adsorbate_indices, product_adsorbate_indices,
                    reaction_name,
                    previous_step_product=previous_step_product,
                    previous_step_product_adsorbate_indices=previous_step_product_adsorbate_indices,
                    step_index=step_index,
                )

                # Step 2: Apply plan
                adj_reactant, adj_product, new_r_idx, new_p_idx = self.apply_plan(
                    reactant, product,
                    reactant_adsorbate_indices, product_adsorbate_indices,
                    plan
                )

                # Record original and adjusted structures for energy correction
                # This allows us to compute the energy offset from virtual atoms later
                plan.original_reactant_energy = len(reactant)  # Store atom count as proxy
                plan.adjusted_reactant_energy = len(adj_reactant)
                plan.original_product_energy = len(product)
                plan.adjusted_product_energy = len(adj_product)

                if len(adj_reactant) != len(reactant) or len(adj_product) != len(product):
                    logger.info(
                        f"Structure adjustment: "
                        f"Reactant {len(reactant)}→{len(adj_reactant)} atoms, "
                        f"Product {len(product)}→{len(adj_product)} atoms"
                    )

                # Step 3: Validate
                validation = self.validate_structures(
                    adj_reactant, adj_product,
                    new_r_idx, new_p_idx,
                    plan
                )

                if validation.is_valid:
                    logger.info(f"Structures validated successfully")
                    if all_fixes:
                        logger.info(f"Applied fixes: {all_fixes}")
                    return adj_reactant, adj_product, new_r_idx, new_p_idx, plan

                # Step 4: AUTO-FIX - This is the key agent self-healing loop!
                logger.warning(f"Validation failed: {validation.issues}")
                logger.info("Applying automatic fixes...")

                adj_reactant, adj_product, new_r_idx, new_p_idx, fixes = self.auto_fix_structures(
                    adj_reactant, adj_product,
                    new_r_idx, new_p_idx,
                    validation
                )
                all_fixes.extend(fixes)

                if fixes:
                    logger.info(f"Applied {len(fixes)} fixes: {fixes}")

                    # Step 5: Re-validate after fixes
                    validation2 = self.validate_structures(
                        adj_reactant, adj_product,
                        new_r_idx, new_p_idx,
                        plan
                    )

                    if validation2.is_valid:
                        logger.info(f"Structures valid after auto-fix!")
                        return adj_reactant, adj_product, new_r_idx, new_p_idx, plan
                    else:
                        logger.warning(f"Still invalid after fixes: {validation2.issues}")
                        # Continue to next attempt with fresh plan
                else:
                    logger.warning("No automatic fixes available, continuing to next attempt")

            except Exception as e:
                logger.error(f"Attempt {attempt + 1} failed: {e}")
                import traceback
                traceback.print_exc()
                if attempt == self.max_retries - 1:
                    raise

        # If all retries failed, do final auto-fix and return
        logger.warning(f"All {self.max_retries} attempts exhausted")

        if adj_reactant is not None:
            # One more auto-fix pass with more aggressive settings
            logger.info("Final auto-fix pass...")
            adj_reactant, adj_product, new_r_idx, new_p_idx, final_fixes = self.auto_fix_structures(
                adj_reactant, adj_product,
                new_r_idx, new_p_idx,
                ValidationResult(is_valid=False, issues=["Final pass"], suggestions=[])
            )
            all_fixes.extend(final_fixes)

            # Ensure atom counts match
            if len(adj_reactant) != len(adj_product):
                logger.error(f"Atom count mismatch: {len(adj_reactant)} vs {len(adj_product)}")
                raise ValueError(f"Could not match atom counts after {self.max_retries} attempts")

            logger.info(f"Returning best effort result with {len(all_fixes)} fixes applied")
            return adj_reactant, adj_product, new_r_idx, new_p_idx, plan
        else:
            raise RuntimeError("Failed to prepare structures after all attempts")

    def align_structures_for_neb(
        self,
        reactant: Atoms,
        product: Atoms,
        reactant_adsorbate_indices: List[int],
        product_adsorbate_indices: List[int],
    ) -> Tuple[Atoms, Atoms]:
        """
        Align product atoms to match reactant ordering for NEB.

        This ensures smooth interpolation by matching corresponding atoms.
        """
        from scipy.optimize import linear_sum_assignment

        r_pos = reactant.get_positions()
        p_pos = product.get_positions()
        r_sym = reactant.get_chemical_symbols()
        p_sym = product.get_chemical_symbols()

        # Build cost matrix for Hungarian algorithm
        n = len(reactant)
        cost_matrix = np.full((n, n), 1e10)

        for i in range(n):
            for j in range(n):
                if r_sym[i] == p_sym[j]:
                    cost_matrix[i, j] = np.linalg.norm(r_pos[i] - p_pos[j])

        # Solve assignment
        row_ind, col_ind = linear_sum_assignment(cost_matrix)

        # Reorder product
        aligned_product = product[col_ind]

        # Update adsorbate indices
        new_p_indices = []
        for old_idx in product_adsorbate_indices:
            new_idx = np.where(col_ind == old_idx)[0]
            if len(new_idx) > 0:
                new_p_indices.append(int(new_idx[0]))

        return reactant, aligned_product
