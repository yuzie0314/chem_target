"""Environment check + ProLIF reference-IFP library builder (3D fallback, Option 1).

Builds ``db/prolif_reference_ifp.json`` (the reference contract consumed by
``utils.fallback_3d.ProLIFFallback``) from PDB co-crystal complexes of the target
class.  Using real crystal poses for the *reference* side means no docking is
needed to build the library — docking is only required at query time
(``ProLIFFallback._dock``).

Usage
-----
    # 1. Check what's installed (safe, no heavy work)
    python utils/build_prolif_reference.py check

    # 2. Build the serine-protease reference IFP library
    python utils/build_prolif_reference.py build --target "serine protease"

Dependencies
------------
Build-time : RDKit, ProLIF, MDAnalysis, requests   (real crystal poses → IFP)
Query-time : a docking backend (smina / gnina / vina) — checked but NOT needed
             to build the reference library.

All heavy imports are lazy so ``check`` runs even when deps are missing.
This is a sanctioned db-writer (like interaction_analyzer.py / pose_extractor.py);
its output db/prolif_reference_ifp.json is gitignored — rebuild offline.

KNOWN ISSUE / TODO (protein prep) — 2026-06-16
----------------------------------------------
The pipeline now: fixes conda OpenBabel's plugin path (``_setup_openbabel``),
protonates the complex, and builds the ProLIF protein from the protonated PDB +
the ligand from a standalone-protonated copy (OpenBabel drops HETATM resnames on
a full-complex rewrite, so the two are sourced separately).  BUT OpenBabel mangles
protein chain/residue topology when it rewrites the protein, so ProLIF then raises
``KeyError: ResidueId(...)`` in ``Fingerprint.generate`` (an interaction references
a protein residue missing from the rebuilt residue table).  **Fix needed: prepare
the protein with a topology-preserving tool — PDBFixer (openmm) or ``reduce`` —
instead of OpenBabel.**  Use OpenBabel only for the small-molecule ligand.
Suggested: ``conda install -c conda-forge pdbfixer openmm`` then replace the
protein branch of ``_ifp_from_complex`` with a PDBFixer addMissingHydrogens pass.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
from pathlib import Path

# Force UTF-8 console so the report's unicode glyphs don't mangle on cp950 (Windows).
for _stream in (sys.stdout, sys.stderr):
    if hasattr(_stream, "reconfigure"):
        _stream.reconfigure(encoding="utf-8")

_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT))   # so `import utils.fallback_3d` works when run as a script
_OUT_PATH = _ROOT / "db" / "prolif_reference_ifp.json"
_PDB_DIR = _ROOT / "db" / "prolif_pdb"          # cached downloaded PDBs (gitignored)
_RCSB_PDB_URL = "https://files.rcsb.org/download/{pdb}.pdb"

# ── Seed reference complexes (PDB co-crystals) per target class ──────────────────
# Starter set — EXPAND with more non-benzamidine peptidomimetic co-crystals to
# cover the 8 serine-protease misses that lack a Benzamidine FG.
# ligand_resname = the het residue name of the bound inhibitor in that PDB.
_SEED_REFERENCES: dict[str, dict] = {
    "serine protease": {
        "binding_site_residues": [
            "HIS57", "ASP102", "SER195",     # catalytic triad
            "ASP189", "SER190", "GLY216", "GLY219",  # S1 pocket
        ],
        "complexes": [
            {"pdb": "3PTB", "ligand_resname": "BEN", "note": "trypsin + benzamidine"},
            {"pdb": "1OYT", "ligand_resname": "FSN", "note": "thrombin + non-covalent inhibitor"},
            {"pdb": "1DWD", "ligand_resname": "MID", "note": "thrombin + NAPAP"},
            {"pdb": "2ZFF", "ligand_resname": "53U", "note": "factor Xa + peptidomimetic"},
            {"pdb": "1F0R", "ligand_resname": "815", "note": "factor Xa + non-benzamidine"},
        ],
    },
}


# ── Environment check ───────────────────────────────────────────────────────────

def _probe_import(modname: str) -> tuple[bool, str]:
    """Return (available, version-or-error) for an importable module."""
    try:
        mod = __import__(modname)
        return True, getattr(mod, "__version__", "ok")
    except Exception as exc:  # noqa: BLE001
        return False, f"{type(exc).__name__}: {exc}"


def check_environment() -> dict:
    """Print a dependency status table; return a dict of results."""
    build_deps = {
        "rdkit":       "core cheminformatics (have it)",
        "prolif":      "interaction fingerprints  →  pip install prolif",
        "MDAnalysis":  "PDB/trajectory parsing     →  pip install MDAnalysis",
        "requests":    "PDB download               →  pip install requests",
    }
    backends = ["smina", "gnina", "vina", "qvina2"]

    print("=" * 68)
    print("  ProLIF 3D-fallback — environment check")
    print("=" * 68)
    print("\nBUILD-TIME (needed to build the reference IFP library):")
    results: dict = {"build": {}, "docking": {}}
    all_build_ok = True
    for mod, hint in build_deps.items():
        ok, info = _probe_import(mod)
        all_build_ok &= ok
        results["build"][mod] = ok
        mark = "OK " if ok else "MISSING"
        print(f"  [{mark:7}] {mod:12} {info if ok else hint}")

    print("\nQUERY-TIME (docking backend — needed later, NOT for building):")
    any_backend = False
    for b in backends:
        path = shutil.which(b)
        results["docking"][b] = bool(path)
        any_backend |= bool(path)
        if path:
            print(f"  [OK     ] {b:12} {path}")
    if not any_backend:
        print("  [MISSING] no docking backend on PATH")
        print("            install one for query-time docking, e.g. "
              "`conda install -c conda-forge smina` (gnina needs GPU/Docker)")

    results["build_ready"] = all_build_ok
    results["docking_ready"] = any_backend
    print("\n" + "-" * 68)
    print(f"  Build reference library: {'READY' if all_build_ok else 'BLOCKED (install build deps)'}")
    print(f"  Query-time docking:      {'READY' if any_backend else 'not yet (ok for now)'}")
    print("=" * 68)
    return results


# ── PDB fetch ────────────────────────────────────────────────────────────────────

def _fetch_pdb(pdb_id: str) -> Path | None:
    """Download a PDB file to the cache; return path or None on failure."""
    _PDB_DIR.mkdir(parents=True, exist_ok=True)
    dest = _PDB_DIR / f"{pdb_id}.pdb"
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    import requests  # lazy
    try:
        r = requests.get(_RCSB_PDB_URL.format(pdb=pdb_id), timeout=30)
        if r.status_code == 200 and r.text.strip():
            dest.write_text(r.text, encoding="utf-8")
            return dest
        print(f"    [warn] {pdb_id}: HTTP {r.status_code}")
    except Exception as exc:  # noqa: BLE001
        print(f"    [warn] {pdb_id}: {exc}")
    return None


# ── IFP from a crystal complex ──────────────────────────────────────────────────

def _setup_openbabel() -> None:
    """Make conda's OpenBabel find its format plugins (Windows needs this).

    The .obf plugins live in <env>/Library/bin and load only when that dir is on
    PATH + the DLL search path, with BABEL_DATADIR/BABEL_LIBDIR set. Paths are
    derived from sys.prefix so this is portable across machines.
    """
    prefix = Path(sys.prefix)
    libbin = prefix / "Library" / "bin"
    data = prefix / "share" / "openbabel"
    if libbin.exists():
        os.environ.setdefault("BABEL_LIBDIR", str(libbin))
        os.environ["PATH"] = str(libbin) + os.pathsep + os.environ.get("PATH", "")
        try:
            os.add_dll_directory(str(libbin))   # Windows DLL resolution
        except (AttributeError, OSError):
            pass
    if data.exists():
        os.environ.setdefault("BABEL_DATADIR", str(data))


def _protonate(pdb_path: Path) -> Path:
    """Add explicit hydrogens to a crystal PDB (ProLIF/RDKit need them).

    Crystal structures lack H, but MDAnalysis->RDKit bond-order inference requires
    explicit H. Protonate the whole complex at pH 7.4 with OpenBabel (pybel),
    cached as {pdb}_H.pdb. Returns the protonated path (or the original on failure).
    """
    out = pdb_path.with_name(pdb_path.stem + "_H.pdb")
    if out.exists() and out.stat().st_size > 0:
        return out
    try:
        _setup_openbabel()
        from openbabel import pybel  # lazy
        mol = next(pybel.readfile("pdb", str(pdb_path)))
        mol.OBMol.AddHydrogens(False, True, 7.4)   # (polar_only, correct_for_pH, pH)
        mol.write("pdb", str(out), overwrite=True)
        return out
    except Exception as exc:  # noqa: BLE001
        print(f"    [warn] protonation failed for {pdb_path.name}: {exc}")
        return pdb_path


def _prep_protein_pdbfixer(pdb_path: Path) -> Optional[Path]:
    """Protonate the protein with PDBFixer (topology-preserving) → cached PDB.

    Unlike OpenBabel (which mangles chain/residue numbering and triggers a ProLIF
    ResidueId KeyError), PDBFixer keeps the protein topology intact: it strips
    heterogens/water, fills missing atoms, and adds hydrogens at pH 7.0. Returns
    the protonated protein PDB path, or None on failure.
    """
    out = pdb_path.with_name(pdb_path.stem + "_pf.pdb")
    if out.exists() and out.stat().st_size > 0:
        return out
    try:
        from pdbfixer import PDBFixer            # lazy
        from openmm.app import PDBFile            # lazy
        fixer = PDBFixer(filename=str(pdb_path))
        fixer.removeHeterogens(keepWater=False)   # protein only; ligand handled separately
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.0)
        with open(out, "w") as fh:
            PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)
        return out
    except Exception as exc:  # noqa: BLE001
        print(f"    [warn] PDBFixer protein prep failed for {pdb_path.name}: {exc}")
        return None


def _protonate_ligand(original_pdb: Path, resname: str):
    """Extract a ligand by resname from the ORIGINAL PDB, protonate it standalone
    (OpenBabel perceives bond orders from the crystal geometry + adds H), and
    return an RDKit mol. Done separately because OpenBabel drops HETATM resnames
    when it rewrites the whole complex.
    """
    import tempfile
    import MDAnalysis as mda          # lazy
    from rdkit import Chem            # lazy
    u = mda.Universe(str(original_pdb))
    sel = u.select_atoms(f"resname {resname}")
    if sel.n_atoms == 0:
        return None
    _setup_openbabel()
    from openbabel import pybel       # lazy
    with tempfile.TemporaryDirectory() as td:
        raw, sdf = Path(td) / "lig.pdb", Path(td) / "lig.sdf"
        sel.write(str(raw))
        m = next(pybel.readfile("pdb", str(raw)))
        m.OBMol.AddHydrogens(False, True, 7.4)
        m.write("sdf", str(sdf), overwrite=True)
        mols = [x for x in Chem.SDMolSupplier(str(sdf), removeHs=False) if x is not None]
        return mols[0] if mols else None


def _ifp_from_complex(pdb_path: Path, ligand_resname: str):
    """Compute the ProLIF interaction-key set for one PDB protein-ligand complex.

    Protein comes from the protonated complex (OpenBabel keeps standard AA
    resnames + adds H); the ligand is extracted from the original PDB and
    protonated standalone (OpenBabel loses HETATM resnames on full-complex
    rewrite). Returns (ifp_keys: list[str], n_keys) or (None, 0). Same
    'PROTRES.Interaction' key representation as utils.fallback_3d. Lazy imports.
    """
    import MDAnalysis as mda          # lazy
    import prolif as plf              # lazy
    from utils.fallback_3d import ifp_keys_from_fingerprint

    prot_pdb = _prep_protein_pdbfixer(pdb_path)   # topology-preserving H (avoids ProLIF KeyError)
    if prot_pdb is None:
        return None, 0
    prot_sel = mda.Universe(str(prot_pdb)).select_atoms("protein")
    if prot_sel.n_atoms == 0:
        print(f"    [warn] {pdb_path.name}: protein not found after PDBFixer prep")
        return None, 0
    lig_rdkit = _protonate_ligand(pdb_path, ligand_resname)
    if lig_rdkit is None:
        print(f"    [warn] {pdb_path.name}: ligand '{ligand_resname}' prep failed")
        return None, 0

    prot = plf.Molecule.from_mda(prot_sel)
    lig = plf.Molecule.from_rdkit(lig_rdkit)
    fp = plf.Fingerprint()
    fp.run_from_iterable([lig], prot, progress=False, n_jobs=1)  # serial: avoids parallel ResidueId KeyError
    keys = ifp_keys_from_fingerprint(fp)
    return (keys, len(keys)) if keys else (None, 0)


# ── Build ────────────────────────────────────────────────────────────────────────

def build_reference(target_class: str, out_path: Path = _OUT_PATH) -> None:
    """Build/refresh the reference IFP library for one target class."""
    if target_class not in _SEED_REFERENCES:
        print(f"No seed complexes for '{target_class}'. "
              f"Known: {list(_SEED_REFERENCES)}", file=sys.stderr)
        sys.exit(1)

    env = check_environment()
    if not env["build_ready"]:
        print("\nBuild blocked — install the MISSING build-time deps above first.",
              file=sys.stderr)
        sys.exit(1)

    spec = _SEED_REFERENCES[target_class]
    print(f"\nBuilding reference IFPs for '{target_class}' "
          f"({len(spec['complexes'])} complexes) ...")
    actives = []
    for c in spec["complexes"]:
        pdb_path = _fetch_pdb(c["pdb"])
        if pdb_path is None:
            continue
        ifp_keys, n_keys = _ifp_from_complex(pdb_path, c["ligand_resname"])
        if ifp_keys is None:
            continue
        actives.append({
            "source_pdb": c["pdb"],
            "ligand_resname": c["ligand_resname"],
            "note": c.get("note", ""),
            "ifp_keys": ifp_keys,
            "n_keys": n_keys,
        })
        print(f"    [ok] {c['pdb']} ({c['ligand_resname']}): {n_keys} interaction keys")

    # Docking receptor for query-time: use the first complex (real binding site).
    first = spec["complexes"][0]

    # Merge into existing library (preserve other target classes)
    library = {}
    if out_path.exists():
        library = json.loads(out_path.read_text(encoding="utf-8"))
    library[target_class] = {
        "binding_site_residues": spec["binding_site_residues"],
        "docking": {"receptor_pdb": first["pdb"],
                    "autobox_ligand_resname": first["ligand_resname"]},
        "actives": actives,
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(library, ensure_ascii=False, indent=2),
                        encoding="utf-8")
    print(f"\nSaved {len(actives)} reference IFPs for '{target_class}' → {out_path}")
    if not actives:
        print("  WARNING: 0 IFPs built — check ligand resnames / network.")


# ── CLI ──────────────────────────────────────────────────────────────────────────

def main() -> None:
    p = argparse.ArgumentParser(description="ProLIF 3D-fallback env check + reference builder")
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("check", help="environment dependency check only")
    b = sub.add_parser("build", help="build the reference IFP library")
    b.add_argument("--target", default="serine protease",
                   help="target class to build (default: serine protease)")
    b.add_argument("--out", default=str(_OUT_PATH), help="output JSON path")
    args = p.parse_args()

    if args.cmd == "check":
        check_environment()
    elif args.cmd == "build":
        build_reference(args.target, Path(args.out))


if __name__ == "__main__":
    main()
