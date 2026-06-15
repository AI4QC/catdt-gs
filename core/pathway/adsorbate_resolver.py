"""Resolve adsorbate labels to ASE Atoms, with fairchem DB + RDKit + LLM fallback.

Lookup order:
1. fairchem adsorbate DB (exact label match)
2. Built-in alias mapping (e.g. *HCOO* -> *OCHO)
3. Built-in SMILES table -> RDKit 3D generation
4. LLM-generated SMILES -> RDKit 3D -> auto-persist to DB
"""

from __future__ import annotations

import json
import logging
import os
import pickle
import re
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from ase import Atoms

logger = logging.getLogger(__name__)

_LABEL_ALIASES: Dict[str, str] = {
    "*CO": "*CO",
    "*H2O": "*OH2",
    "*H2O*": "*OH2",
    "*OCH3*": "*OCH3",
    "*CH3O*": "*OCH3",
    "*HCOO*": "*OCHO",
    "*HCOO": "*OCHO",
    "*COOH*": "*COOH",
    "*H2CO*": "*CH2*O",
    "*H2CO": "*CH2*O",
    "*CH2O*": "*CH2*O",
    "*CHO*": "*CHO",
    "*H3CO*": "*OCH3",
    "*H3CO": "*OCH3",
    "*CH3OH": "*OHCH3",
    "*CH3OH*": "*OHCH3",
    "*CH2OH*": "*CH2OH",
    "*CH2OH": "*CH2OH",
    "*CHOH*": "*CHOH",
    "*OH*": "*OH",
    "*OOH*": "*OOH",
}

_SMILES_TABLE: Dict[str, Tuple[str, List[int]]] = {
    # PDH hydrocarbon chain fallback. The fairchem DB is the source of truth;
    # keep these templates for missing entries and DB repair.
    "*C3H8": ("CCC", [0]),
    "*C3H7": ("[CH2]CC", [0]),
    "*C3H6": ("C=CC", [0]),
    "*C3H5": ("[CH2]C=C", [0]),
    "*C2H5": ("[CH2]C", [0]),
    # CO2: C binds to surface
    "*CO2": ("O=C=O", [1]),
    "*CO2*": ("O=C=O", [1]),
    # HCOOH (formic acid): carbonyl O binds to surface
    "*HCOOH*": ("OC=O", [0]),
    "*HCOOH": ("OC=O", [0]),
    # H2COOH: radical-like CH3O2, O binds to surface
    "*H2COOH*": ("O[CH]O", [0]),
    "*H2COOH": ("O[CH]O", [0]),
    "*HCOH*": ("OC", [0]),
    "*H2CO2*": ("OC(=O)", [0]),
}


def _load_fairchem_db() -> Tuple[Dict[int, tuple], str]:
    from fairchem.data.oc.databases.pkls import ADSORBATE_PKL_PATH
    with open(ADSORBATE_PKL_PATH, "rb") as f:
        db = pickle.load(f)
    return db, ADSORBATE_PKL_PATH


def _db_label_set(db: Dict[int, tuple]) -> Dict[str, int]:
    return {db[k][1]: k for k in db}


def _smiles_to_atoms(smiles: str, binding_indices: Optional[List[int]] = None) -> Tuple[Atoms, List[int]]:
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit cannot parse SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    status = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    if status != 0:
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
    except Exception:
        pass

    conf = mol.GetConformer()
    symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
    positions = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])

    atoms = Atoms(symbols=symbols, positions=positions)

    if binding_indices is None:
        heavy = [i for i, s in enumerate(symbols) if s != "H"]
        binding_indices = [heavy[0]] if heavy else [0]

    return atoms, binding_indices


def _llm_generate_smiles(label: str) -> Tuple[str, List[int]]:
    import sys
    sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[2]))
    from camel_agents.camel_model_backend import get_camel_model_backend

    backend = get_camel_model_backend(temperature=0.0, max_tokens=300)

    clean = re.sub(r"^\*|\*$", "", label).strip()

    prompt = (
        f"I am building an adsorbate molecule for DFT slab calculations. "
        f"The surface-science label is '{label}' (species: '{clean}').\n\n"
        f"This molecule will be placed on a metal catalyst surface. "
        f"The binding atom(s) (marked by * in the label) sit closest to the surface (lowest z). "
        f"The rest of the molecule extends upward away from the surface.\n\n"
        f"Rules:\n"
        f"1. Give the SMILES for the NEUTRAL, CLOSED-SHELL molecule (saturate with H if needed).\n"
        f"2. binding_indices: 0-based heavy-atom indices of the atom(s) that bond to the surface. "
        f"In the label, * at the beginning means the FIRST heavy atom binds; * at the end means the LAST.\n"
        f"3. The molecular formula must match: only C, H, O (no N, S, metals).\n\n"
        f"Examples:\n"
        f'  *OCHO  -> {{"smiles":"O=CO","binding_indices":[2],"formula":"CHO2"}}  (C binds to surface via the terminal O)\n'
        f'  *CO2   -> {{"smiles":"O=C=O","binding_indices":[1],"formula":"CO2"}}   (C binds)\n'
        f'  *OCH3  -> {{"smiles":"CO","binding_indices":[1],"formula":"CH4O"}}     (O binds, CH3 points up)\n'
        f'  *OHCH3 -> {{"smiles":"CO","binding_indices":[1],"formula":"CH4O"}}     (O binds)\n'
        f'  *CH2*O -> {{"smiles":"C=O","binding_indices":[0,1],"formula":"CH2O"}}  (both C and O bind)\n\n'
        f"Return ONLY the JSON object, no explanation:\n"
        f'{{"smiles": "...", "binding_indices": [...], "formula": "..."}}'
    )

    logger.info("Calling LLM to generate SMILES for '%s'...", label)

    response = backend.run([{"role": "user", "content": prompt}])
    content = response.choices[0].message.content or ""
    content = content.strip()

    if not content:
        raise RuntimeError(f"LLM returned empty content for '{label}'")

    content = re.sub(r"^```(?:json)?\s*", "", content)
    content = re.sub(r"\s*```$", "", content)

    data = json.loads(content)
    smiles = str(data["smiles"])
    binding = [int(x) for x in data.get("binding_indices", [0])]

    from rdkit import Chem
    if Chem.MolFromSmiles(smiles) is None:
        raise ValueError(f"LLM returned invalid SMILES: {smiles}")

    logger.info("LLM generated SMILES for '%s': %s (binding=%s)", label, smiles, binding)
    return smiles, binding


def _persist_to_db(
    db: Dict[int, tuple],
    db_path: str,
    atoms: Atoms,
    label: str,
    binding_indices: List[int],
) -> int:
    new_id = max(db.keys()) + 1
    rxn_str = f"* &rarr; {label}"
    db[new_id] = (atoms.copy(), label, np.array(binding_indices), rxn_str)
    with open(db_path, "wb") as f:
        pickle.dump(db, f)
    logger.info("Persisted adsorbate '%s' to DB as id=%d (path=%s)", label, new_id, db_path)
    return new_id


def resolve_adsorbate(label: str) -> Tuple[Atoms, List[int]]:
    """Resolve an adsorbate label to (Atoms, binding_indices).

    Tries fairchem DB -> alias -> SMILES table -> LLM, persisting new
    entries to the DB for future runs.
    """
    db, db_path = _load_fairchem_db()
    label_map = _db_label_set(db)

    if label in label_map:
        entry = db[label_map[label]]
        return entry[0].copy(), list(entry[2])

    alias = _LABEL_ALIASES.get(label)
    if alias and alias in label_map:
        logger.info("Resolved '%s' via alias -> '%s'", label, alias)
        entry = db[label_map[alias]]
        return entry[0].copy(), list(entry[2])

    if label in _SMILES_TABLE:
        smiles, bind_idx = _SMILES_TABLE[label]
        logger.info("Building '%s' from built-in SMILES: %s", label, smiles)
        atoms, binding = _smiles_to_atoms(smiles, bind_idx)
        _persist_to_db(db, db_path, atoms, label, binding)
        return atoms, binding

    logger.warning("Adsorbate '%s' not in DB/aliases/SMILES table — calling LLM", label)
    smiles, bind_idx = _llm_generate_smiles(label)
    atoms, binding = _smiles_to_atoms(smiles, bind_idx)
    _persist_to_db(db, db_path, atoms, label, binding)
    return atoms, binding
