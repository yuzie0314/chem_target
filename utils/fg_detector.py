"""Functional group detection logic using RDKit SMARTS patterns.

Primary API
-----------
detect_smarts(smiles)          → list[str]          single molecule, FG presence
detect_smarts_table(compounds) → pd.DataFrame       multi-compound abundance table

Both use the same FG_SMARTS patterns as interaction_analyzer.py, plus the
Python-based detector for Steroid scaffold, keeping the FG profile consistent
end-to-end (query compound ↔ BioLiP residue table).

Steroid detection
-----------------
The steroid scaffold (6-6-6-5 fused tetracyclic) cannot be expressed as a
valid RDKit SMARTS because the ``rN`` ring-size primitive is computed from the
SSSR, and the ring-junction atoms between the 5- and 6-membered D/C rings are
assigned to the smallest ring only (r5), so ``[r5;r6]`` never matches.
Detection uses ``_detect_steroid_core(mol)`` instead — see docstring below.

Legacy API (kept for backward compatibility, not used in main pipeline)
-----------------------------------------------------------------------
detect(compounds) → pd.DataFrame  RDKit fr_* fragment counters
"""

from __future__ import annotations

from collections import defaultdict
from typing import Callable

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


# ── Steroid scaffold Python-based detector ─────────────────────────────────────

def _detect_steroid_core(mol: Chem.Mol) -> bool:
    """Detect steroid 6-6-6-5 fused tetracyclic scaffold.

    Algorithm
    ---------
    For each *connected ring system* in the molecule, counts:
    - ``r5_C``: carbon atoms in a 5-membered ring
    - ``r6_C``: carbon atoms in a 6-membered ring
    - ``both_C``: carbon atoms in *both* a 5- and a 6-membered ring
                  (the two C-D ring junction atoms in steroid numbering)

    A connected ring system is classified as steroid-like when:
        r5_C >= 5   (a complete 5-membered ring, all carbons)
        r6_C >= 10  (at least two 6-membered rings)
        both_C >= 2 (the two C-D junction atoms present)

    These three criteria collectively distinguish the tetracyclic steroid ABCD
    scaffold from bicyclic (indane: r6_C=6) and purely acyclic structures.
    Vitamin D secosteroids and similar ring-opened analogues are correctly
    excluded (r6_C < 10 or both_C < 2).

    Why not SMARTS?
    ---------------
    RDKit's ``rN`` SMARTS primitive uses the SSSR (Smallest Set of Smallest
    Rings). For steroids whose SSSR happens to capture only 3 rings (e.g. the
    SSSR picks A, B, D but expresses C as a combination), the C-D junction
    atoms are assigned exclusively to the 5-membered ring (r5), making
    ``[r5;r6]`` always fail.  ``IsAtomInRingOfSize`` does not depend on the
    SSSR and reliably identifies ring membership.

    Args:
        mol: RDKit molecule object (caller must ensure mol is not None).

    Returns:
        True if at least one connected ring system matches steroid criteria.
    """
    ri = mol.GetRingInfo()
    all_ring_atoms: set[int] = set()
    for ring in ri.AtomRings():
        all_ring_atoms.update(ring)
    if not all_ring_atoms:
        return False

    # Build adjacency within ring system (only ring bonds)
    adj: dict[int, set[int]] = defaultdict(set)
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in all_ring_atoms and b in all_ring_atoms and bond.IsInRing():
            adj[a].add(b)
            adj[b].add(a)

    # Find connected ring components via BFS
    visited: set[int] = set()
    ring_systems: list[set[int]] = []
    for start in all_ring_atoms:
        if start in visited:
            continue
        component: set[int] = set()
        queue = [start]
        while queue:
            curr = queue.pop()
            if curr in visited:
                continue
            visited.add(curr)
            component.add(curr)
            queue.extend(adj[curr] - visited)
        ring_systems.append(component)

    # Evaluate each connected ring system
    for sys_atoms in ring_systems:
        r5_C = sum(
            1 for i in sys_atoms
            if mol.GetAtomWithIdx(i).GetAtomicNum() == 6
            and ri.IsAtomInRingOfSize(i, 5)
        )
        r6_C = sum(
            1 for i in sys_atoms
            if mol.GetAtomWithIdx(i).GetAtomicNum() == 6
            and ri.IsAtomInRingOfSize(i, 6)
        )
        both_C = sum(
            1 for i in sys_atoms
            if mol.GetAtomWithIdx(i).GetAtomicNum() == 6
            and ri.IsAtomInRingOfSize(i, 5)
            and ri.IsAtomInRingOfSize(i, 6)
        )
        if r5_C >= 5 and r6_C >= 10 and both_C >= 2:
            return True

    return False


# ── Python-based FG detectors (for FGs that cannot use SMARTS) ───────────────

_PYTHON_DETECTORS: dict[str, Callable[[Chem.Mol], bool]] = {
    "Steroid": _detect_steroid_core,
}


# ── SMARTS-based detection (primary) ─────────────────────────────────────────

def detect_smarts(smiles: str) -> list[str]:
    """Detect functional groups in a single SMILES via SMARTS + Python detectors.

    Uses FG_SMARTS patterns for all FGs except "Steroid", which is detected by
    ``_detect_steroid_core``.  Consistent with the FG × residue table built by
    interaction_analyzer.py.

    Args:
        smiles: Input molecule as SMILES string.

    Returns:
        List of FG names (keys of FG_SMARTS plus "Steroid") present in the
        molecule.  Empty list if SMILES is invalid or no FG matches.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []

    results: list[str] = [
        name
        for name, pattern in _SMARTS_PATTERNS.items()
        if pattern is not None and mol.GetSubstructMatches(pattern)
    ]

    # Python-based detectors (e.g. Steroid)
    for fg_name, detector in _PYTHON_DETECTORS.items():
        if detector(mol):
            results.append(fg_name)

    return results


def detect_smarts_table(compounds: dict[str, str]) -> pd.DataFrame:
    """Build a FG abundance table for multiple compounds.

    Counts SMARTS substructure matches (via GetSubstructMatches) plus calls
    Python-based detectors (e.g. Steroid, which returns 0 or 1).

    Rows with all-zero counts are dropped so the table only shows FGs that
    appear in at least one compound.

    Args:
        compounds: {compound_name: smiles}

    Returns:
        DataFrame with fg_name index, one column per compound, integer counts.
        Rows with all-zero counts are dropped.
    """
    rows: dict[str, dict[str, int]] = {}

    # SMARTS-based FGs
    for fg_name, pattern in _SMARTS_PATTERNS.items():
        row: dict[str, int] = {}
        for comp_name, smiles in compounds.items():
            mol = Chem.MolFromSmiles(smiles)
            row[comp_name] = len(mol.GetSubstructMatches(pattern)) if mol else 0
        rows[fg_name] = row

    # Python-based FGs (e.g. Steroid) — count is 0 or 1
    for fg_name, detector in _PYTHON_DETECTORS.items():
        row = {}
        for comp_name, smiles in compounds.items():
            mol = Chem.MolFromSmiles(smiles)
            row[comp_name] = int(detector(mol)) if mol else 0
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
