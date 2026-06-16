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
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

# Force UTF-8 console so the report's unicode glyphs don't mangle on cp950 (Windows).
for _stream in (sys.stdout, sys.stderr):
    if hasattr(_stream, "reconfigure"):
        _stream.reconfigure(encoding="utf-8")

_ROOT = Path(__file__).resolve().parent.parent
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
            {"pdb": "1DWD", "ligand_resname": "MIT", "note": "thrombin + NAPAP"},
            {"pdb": "2ZFF", "ligand_resname": "DX9", "note": "factor Xa + peptidomimetic"},
            {"pdb": "1F0R", "ligand_resname": "RPR", "note": "factor Xa + non-benzamidine"},
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

def _ifp_from_complex(pdb_path: Path, ligand_resname: str):
    """Compute a ProLIF interaction fingerprint for one PDB protein-ligand complex.

    Returns (ifp_hex, n_bits) or (None, 0) on failure. Lazy-imports the heavy stack.
    """
    import MDAnalysis as mda          # lazy
    import prolif as plf              # lazy

    u = mda.Universe(str(pdb_path))
    lig_sel = u.select_atoms(f"resname {ligand_resname}")
    prot_sel = u.select_atoms("protein")
    if lig_sel.n_atoms == 0 or prot_sel.n_atoms == 0:
        print(f"    [warn] {pdb_path.name}: ligand '{ligand_resname}' or protein not found")
        return None, 0

    lig = plf.Molecule.from_mda(lig_sel)
    prot = plf.Molecule.from_mda(prot_sel)
    fp = plf.Fingerprint()
    fp.run_from_iterable([lig], prot, progress=False)
    bvs = fp.to_bitvectors()          # list[ExplicitBitVect], one per ligand frame
    if not bvs:
        return None, 0
    from rdkit.DataStructs import BitVectToText  # lazy
    bv = bvs[0]
    return bv.ToBase64() if hasattr(bv, "ToBase64") else BitVectToText(bv), bv.GetNumBits()


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
        ifp_hex, n_bits = _ifp_from_complex(pdb_path, c["ligand_resname"])
        if ifp_hex is None:
            continue
        actives.append({
            "source_pdb": c["pdb"],
            "ligand_resname": c["ligand_resname"],
            "note": c.get("note", ""),
            "ifp_hex": ifp_hex,
            "n_bits": n_bits,
        })
        print(f"    [ok] {c['pdb']} ({c['ligand_resname']}): {n_bits}-bit IFP")

    # Merge into existing library (preserve other target classes)
    library = {}
    if out_path.exists():
        library = json.loads(out_path.read_text(encoding="utf-8"))
    library[target_class] = {
        "binding_site_residues": spec["binding_site_residues"],
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
