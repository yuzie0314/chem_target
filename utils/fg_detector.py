"""Functional group detection logic using RDKit fragment counters."""

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Fragments

from constants.fg_names import FG_NAMES

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
