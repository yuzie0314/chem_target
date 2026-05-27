"""Functional group detection logic using RDKit fragment counters and SMARTS patterns."""

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Fragments

from constants.fg_names import FG_NAMES
from constants.fg_smarts import FG_SMARTS

# All callable fr_* functions from RDKit Fragments module
_FG_FUNCTIONS: dict[str, object] = {
    name: func
    for name, func in Fragments.__dict__.items()
    if callable(func) and name.startswith("fr_")
}


def detect(compounds: dict[str, str]) -> pd.DataFrame:
    """Detect functional groups across a set of compounds.

    Args:
        compounds: {compound_name: smiles}

    Returns:
        DataFrame with fg_code index, 'functional_group' name column,
        and one column per compound containing integer counts.
        Rows with all-zero counts are dropped.
    """
    rows: dict[str, dict[str, int]] = {}

    for fg_code, fg_func in _FG_FUNCTIONS.items():
        row: dict[str, int] = {}
        for comp_name, smiles in compounds.items():
            mol = Chem.MolFromSmiles(smiles)
            row[comp_name] = fg_func(mol) if mol else 0
        rows[fg_code] = row

    df = pd.DataFrame(rows).T
    df.index.name = "fg_code"

    # Drop functional groups that appear in no compound
    df = df[(df != 0).any(axis=1)]

    # Insert human-readable name as first column
    df.insert(0, "functional_group", df.index.map(lambda x: FG_NAMES.get(x, x)))

    return df


# ── SMARTS-based detection (used by target_predictor) ─────────────────────────

# Pre-compile SMARTS patterns once at import time
_SMARTS_PATTERNS: dict[str, object] = {
    name: Chem.MolFromSmarts(smarts)
    for name, smarts in FG_SMARTS.items()
    if Chem.MolFromSmarts(smarts) is not None
}


def detect_smarts(smiles: str) -> list[str]:
    """Detect functional groups in a SMILES using SMARTS substructure matching.

    Returns a list of FG names (keys of FG_SMARTS) present in the molecule.
    Uses substructure matching rather than RDKit fragment counters — consistent
    with the FG × residue table built by interaction_analyzer.py.

    Args:
        smiles: Input molecule as SMILES string.

    Returns:
        List of matching FG names from FG_SMARTS. Empty list if SMILES is invalid.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []
    return [
        name
        for name, pattern in _SMARTS_PATTERNS.items()
        if pattern is not None and mol.GetSubstructMatches(pattern)
    ]
