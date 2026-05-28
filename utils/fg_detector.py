"""Functional group detection logic using RDKit SMARTS patterns.

Primary API
-----------
detect_smarts(smiles)          → list[str]          single molecule, FG presence
detect_smarts_table(compounds) → pd.DataFrame       multi-compound abundance table

Both use the same FG_SMARTS patterns as interaction_analyzer.py, keeping
the FG profile consistent end-to-end (query compound ↔ BioLiP residue table).

Legacy API (kept for backward compatibility, not used in main pipeline)
-----------------------------------------------------------------------
detect(compounds) → pd.DataFrame  RDKit fr_* fragment counters
"""

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Fragments

from constants.fg_names import FG_NAMES
from constants.fg_smarts import FG_SMARTS

# Suppress RDKit parse noise — invalid SMILES are handled by explicit None checks.
RDLogger.DisableLog("rdApp.error")
RDLogger.DisableLog("rdApp.warning")

# Pre-compile SMARTS patterns once at import time
_SMARTS_PATTERNS: dict[str, Chem.Mol] = {
    name: Chem.MolFromSmarts(smarts)
    for name, smarts in FG_SMARTS.items()
    if Chem.MolFromSmarts(smarts) is not None
}


# ── SMARTS-based detection (primary) ─────────────────────────────────────────

def detect_smarts(smiles: str) -> list[str]:
    """Detect functional groups in a single SMILES via SMARTS substructure matching.

    Consistent with the FG × residue table built by interaction_analyzer.py —
    uses the same FG_SMARTS patterns.

    Args:
        smiles: Input molecule as SMILES string.

    Returns:
        List of FG names (keys of FG_SMARTS) present in the molecule.
        Empty list if SMILES is invalid or no FG matches.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []
    return [
        name
        for name, pattern in _SMARTS_PATTERNS.items()
        if pattern is not None and mol.GetSubstructMatches(pattern)
    ]


def detect_smarts_table(compounds: dict[str, str]) -> pd.DataFrame:
    """Build a FG abundance table for multiple compounds using SMARTS matching.

    Each cell = number of non-overlapping substructure match sites in that
    molecule (via GetSubstructMatches). Rows with all-zero counts are dropped.
    Consistent with the FG × residue table built by interaction_analyzer.py.

    Args:
        compounds: {compound_name: smiles}

    Returns:
        DataFrame with fg_name index, one column per compound, integer counts.
        Rows with all-zero counts are dropped.
    """
    rows: dict[str, dict[str, int]] = {}

    for fg_name, pattern in _SMARTS_PATTERNS.items():
        row: dict[str, int] = {}
        for comp_name, smiles in compounds.items():
            mol = Chem.MolFromSmiles(smiles)
            row[comp_name] = len(mol.GetSubstructMatches(pattern)) if mol else 0
        rows[fg_name] = row

    df = pd.DataFrame(rows).T
    df.index.name = "fg_name"

    # Drop functional groups that appear in no compound
    df = df[(df != 0).any(axis=1)]

    return df


# ── Legacy API: RDKit fr_* fragment counters ──────────────────────────────────

# All callable fr_* functions from RDKit Fragments module
_FG_FUNCTIONS: dict[str, object] = {
    name: func
    for name, func in Fragments.__dict__.items()
    if callable(func) and name.startswith("fr_")
}


def detect(compounds: dict[str, str]) -> pd.DataFrame:
    """[Legacy] Detect functional groups using RDKit fr_* fragment counters.

    Kept for backward compatibility. New code should use detect_smarts_table()
    which uses FG_SMARTS and is consistent with the BioLiP residue table.

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
