"""Input file reading utilities.

Supported formats:
    CSV  — col 0 = compound name, col 1 = SMILES
    SDF  — molecule title as name; SMILES converted via OpenBabel

Planned (extend read_file / add format-specific readers):
    MOL2, InChI strings, SMILES files
"""

import os
import pandas as pd
from rdkit import Chem


# ── CSV ───────────────────────────────────────────────────────────────────────

def read_csv(filepath: str) -> dict[str, str]:
    """Read a CSV file and return a dict of {compound_name: smiles}.

    Expects col 0 = compound name, col 1 = SMILES. Column headers are flexible.
    """
    df = pd.read_csv(filepath)
    df.columns = [c.strip() for c in df.columns]
    name_col = df.columns[0]
    smiles_col = df.columns[1]
    return dict(zip(df[name_col], df[smiles_col].astype(str)))


# ── SDF ───────────────────────────────────────────────────────────────────────

def read_sdf(
    filepath: str,
    name_property: str | None = None,
) -> dict[str, str]:
    """Read an SDF file and return a dict of {compound_name: smiles}.

    Compound names are taken from (in priority order):
      1. The SDF property given by `name_property` (if provided and present)
      2. The molecule title line ($$$$-block title)
      3. A fallback index label: "mol_1", "mol_2", ...

    SMILES are generated from the 2D/3D coordinates via OpenBabel.
    Molecules that OpenBabel cannot parse are skipped with a warning.

    Args:
        filepath:      Path to the .sdf file.
        name_property: Optional SDF property tag to use as the compound name
                       (e.g. "ChEMBL_ID", "CAS", "ID"). Case-sensitive.

    Returns:
        {compound_name: smiles}
    """
    try:
        from openbabel import pybel  # noqa: PLC0415
    except ImportError:
        raise ImportError(
            "OpenBabel (pybel) is required for SDF reading. "
            "Install with: conda install -c conda-forge openbabel"
        )

    compounds: dict[str, str] = {}
    seen_names: dict[str, int] = {}   # track duplicates

    for idx, mol in enumerate(pybel.readfile("sdf", filepath), start=1):
        # ── Name resolution ────────────────────────────────────────────────
        name: str | None = None

        if name_property:
            name = mol.data.get(name_property, "").strip() or None

        if not name:
            title = mol.title.strip()
            name = title if title else None

        if not name:
            name = f"mol_{idx}"

        # De-duplicate: append suffix if name already seen
        if name in seen_names:
            seen_names[name] += 1
            name = f"{name}_{seen_names[name]}"
        else:
            seen_names[name] = 1

        # ── SMILES conversion via OpenBabel ────────────────────────────────
        raw_smi = mol.write("smi").strip()
        if not raw_smi:
            print(f"  ! [{name}] OpenBabel produced empty SMILES — skipped")
            continue

        # Take only the SMILES token (OpenBabel appends the title after a tab)
        smiles = raw_smi.split()[0]

        # ── Validate with RDKit ────────────────────────────────────────────
        if Chem.MolFromSmiles(smiles) is None:
            print(f"  ! [{name}] RDKit could not parse SMILES '{smiles[:60]}' — skipped")
            continue

        compounds[name] = smiles

    return compounds


# ── Auto-dispatch ─────────────────────────────────────────────────────────────

def read_file(
    filepath: str,
    name_property: str | None = None,
) -> dict[str, str]:
    """Read a compound file and return {compound_name: smiles}.

    Dispatches to the correct reader based on file extension:
      .csv         → read_csv
      .sdf / .sd   → read_sdf
      other        → raises ValueError

    Args:
        filepath:      Path to the input file.
        name_property: Passed to read_sdf when reading SDF files (ignored for CSV).
    """
    ext = os.path.splitext(filepath)[1].lower()

    if ext == ".csv":
        return read_csv(filepath)

    if ext in {".sdf", ".sd"}:
        return read_sdf(filepath, name_property=name_property)

    raise ValueError(
        f"Unsupported file format: '{ext}'. "
        "Supported: .csv, .sdf, .sd"
    )
