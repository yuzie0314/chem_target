"""Input file reading utilities → {compound_name: smiles}.

Supported formats:
    CSV          — col 0 = compound name, col 1 = SMILES
    SDF          — RDKit SDMolSupplier; name from property/title; SMILES via RDKit
    MOL2         — OpenBabel (pybel); name from title; SMILES via OpenBabel
    SMI/SMILES   — whitespace-delimited "<SMILES> [name]" per line
    InChI        — "InChI=... [name]" per line; converted via RDKit MolFromInchi

Format choice
-------------
SDF is read with **RDKit** (native, robust, and consistent with the rest of the
RDKit-based pipeline) rather than OpenBabel — conda's OpenBabel format plugins do
not load on Windows without a DLL/PATH fix (see ``_setup_openbabel``), so RDKit is
the reliable path for the common SDF case.  OpenBabel is reserved for MOL2, where
RDKit's reader is weaker; ``_setup_openbabel`` makes its plugins discoverable.

"""

import os
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")   # silence RDKit parse noise; handled via None checks


# ── CSV ───────────────────────────────────────────────────────────────────────

def read_csv(filepath: str) -> dict[str, str]:
    """Read a CSV file and return {compound_name: smiles}.

    Expects col 0 = compound name, col 1 = SMILES. Column headers are flexible.
    """
    df = pd.read_csv(filepath)
    df.columns = [c.strip() for c in df.columns]
    name_col, smiles_col = df.columns[0], df.columns[1]
    return dict(zip(df[name_col], df[smiles_col].astype(str)))


# ── Name de-duplication helper ──────────────────────────────────────────────────

def _dedup(name: str, seen: dict[str, int]) -> str:
    """Return a unique name, appending _2/_3… when a name repeats."""
    if name in seen:
        seen[name] += 1
        return f"{name}_{seen[name]}"
    seen[name] = 1
    return name


# ── SDF (RDKit) ─────────────────────────────────────────────────────────────────

def read_sdf(filepath: str, name_property: str | None = None) -> dict[str, str]:
    """Read an SDF file with RDKit and return {compound_name: smiles}.

    Name priority: (1) the SDF property ``name_property`` if given and present,
    (2) the molecule title (``_Name``), (3) a fallback ``mol_<i>``. Duplicates get
    a numeric suffix. Unparseable records are skipped with a warning.

    Args:
        filepath:      Path to the .sdf / .sd file.
        name_property: Optional SDF property tag to use as the name (e.g. "ChEMBL_ID").
    """
    compounds: dict[str, str] = {}
    seen: dict[str, int] = {}

    supplier = Chem.SDMolSupplier(filepath, removeHs=True, sanitize=True)
    for idx, mol in enumerate(supplier, start=1):
        if mol is None:
            print(f"  ! [mol_{idx}] RDKit could not parse SDF record — skipped")
            continue

        name = None
        if name_property and mol.HasProp(name_property):
            name = mol.GetProp(name_property).strip() or None
        if not name and mol.HasProp("_Name"):
            name = mol.GetProp("_Name").strip() or None
        if not name:
            name = f"mol_{idx}"
        name = _dedup(name, seen)

        smiles = Chem.MolToSmiles(mol)
        if not smiles:
            print(f"  ! [{name}] empty SMILES from RDKit — skipped")
            continue
        compounds[name] = smiles

    return compounds


# ── OpenBabel plugin setup (Windows conda) ──────────────────────────────────────

def _setup_openbabel() -> None:
    """Make conda's OpenBabel find its format plugins (Windows needs this).

    The .obf plugins live in <env>/Library/bin and load only when that dir is on
    PATH + the DLL search path, with BABEL_LIBDIR/BABEL_DATADIR set. Derived from
    sys.prefix → portable across machines.
    """
    prefix = Path(sys.prefix)
    libbin, data = prefix / "Library" / "bin", prefix / "share" / "openbabel"
    if libbin.exists():
        os.environ.setdefault("BABEL_LIBDIR", str(libbin))
        os.environ["PATH"] = str(libbin) + os.pathsep + os.environ.get("PATH", "")
        try:
            os.add_dll_directory(str(libbin))
        except (AttributeError, OSError):
            pass
    if data.exists():
        os.environ.setdefault("BABEL_DATADIR", str(data))


# ── MOL2 (OpenBabel) ─────────────────────────────────────────────────────────────

def read_mol2(filepath: str) -> dict[str, str]:
    """Read a MOL2 file via OpenBabel and return {compound_name: smiles}.

    RDKit's MOL2 reader is fragile; OpenBabel handles MOL2 robustly. Name from the
    molecule title, else ``mol_<i>``. SMILES validated with RDKit.
    """
    _setup_openbabel()
    try:
        from openbabel import pybel  # noqa: PLC0415
    except ImportError as exc:
        raise ImportError(
            "OpenBabel (pybel) is required for MOL2 reading. "
            "Install with: conda install -c conda-forge openbabel"
        ) from exc

    compounds: dict[str, str] = {}
    seen: dict[str, int] = {}
    for idx, mol in enumerate(pybel.readfile("mol2", filepath), start=1):
        name = (mol.title or "").strip() or f"mol_{idx}"
        name = _dedup(name, seen)
        raw = mol.write("smi").strip()
        if not raw:
            print(f"  ! [{name}] OpenBabel produced empty SMILES — skipped")
            continue
        smiles = raw.split()[0]
        if Chem.MolFromSmiles(smiles) is None:
            print(f"  ! [{name}] RDKit could not parse SMILES '{smiles[:60]}' — skipped")
            continue
        compounds[name] = smiles
    return compounds


# ── Plain SMILES file (.smi / .smiles) ──────────────────────────────────────────

def read_smiles_file(filepath: str) -> dict[str, str]:
    """Read a whitespace-delimited SMILES file → {compound_name: smiles}.

    Each non-blank line: ``<SMILES> [name…]`` (the standard .smi convention).
    A leading header line (first token "smiles"/"smi") is skipped. Lines whose
    first token is not a valid SMILES are skipped with a warning. Missing names
    fall back to ``mol_<i>``; duplicates get a numeric suffix.
    """
    compounds: dict[str, str] = {}
    seen: dict[str, int] = {}
    with open(filepath, encoding="utf-8") as fh:
        for idx, line in enumerate(fh, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            tokens = line.split()
            smi = tokens[0]
            if idx == 1 and smi.lower() in {"smiles", "smi"}:
                continue   # header row
            if Chem.MolFromSmiles(smi) is None:
                print(f"  ! [line {idx}] invalid SMILES '{smi[:60]}' — skipped")
                continue
            name = " ".join(tokens[1:]).strip() or f"mol_{idx}"
            compounds[_dedup(name, seen)] = smi
    return compounds


# ── InChI file (.inchi) ──────────────────────────────────────────────────────────

def read_inchi_file(filepath: str) -> dict[str, str]:
    """Read a file of InChI strings → {compound_name: smiles}.

    Each non-blank line: ``InChI=... [name…]`` (InChI strings contain no spaces,
    so the first whitespace token is the InChI, the rest is the name). Converted
    to canonical SMILES via RDKit. Missing names fall back to ``mol_<i>``.
    """
    compounds: dict[str, str] = {}
    seen: dict[str, int] = {}
    with open(filepath, encoding="utf-8") as fh:
        for idx, line in enumerate(fh, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            tokens = line.split()
            inchi = tokens[0]
            if not inchi.startswith("InChI="):
                print(f"  ! [line {idx}] not an InChI string — skipped")
                continue
            mol = Chem.MolFromInchi(inchi)
            if mol is None:
                print(f"  ! [line {idx}] RDKit could not parse InChI — skipped")
                continue
            name = " ".join(tokens[1:]).strip() or f"mol_{idx}"
            compounds[_dedup(name, seen)] = Chem.MolToSmiles(mol)
    return compounds


# ── Auto-dispatch ─────────────────────────────────────────────────────────────

def read_file(filepath: str, name_property: str | None = None) -> dict[str, str]:
    """Read a compound file → {compound_name: smiles}, dispatched by extension.

    .csv → read_csv · .sdf/.sd → read_sdf · .mol2/.mol → read_mol2 · else ValueError.

    Args:
        filepath:      Path to the input file.
        name_property: Passed to read_sdf for SDF files (ignored otherwise).
    """
    ext = os.path.splitext(filepath)[1].lower()
    if ext == ".csv":
        return read_csv(filepath)
    if ext in {".sdf", ".sd"}:
        return read_sdf(filepath, name_property=name_property)
    if ext in {".mol2", ".mol"}:
        return read_mol2(filepath)
    if ext in {".smi", ".smiles"}:
        return read_smiles_file(filepath)
    if ext in {".inchi", ".ich"}:
        return read_inchi_file(filepath)
    raise ValueError(
        f"Unsupported file format: '{ext}'. "
        "Supported: .csv, .sdf, .sd, .mol2, .mol, .smi, .smiles, .inchi"
    )
