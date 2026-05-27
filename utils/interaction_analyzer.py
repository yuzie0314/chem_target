"""Build a functional group × protein residue interaction frequency table.

Data source: BioLiP 2.0 non-redundant flat file (BioLiP_nr.txt.gz)
  https://zhanggroup.org/BioLiP/download/BioLiP_nr.txt.gz

Pipeline:
  1. Download BioLiP_nr (if not cached locally)
  2. Parse: (ligand_3letter_code, binding_residues_1letter[])
  3. Convert ligand 3-letter → SMILES via RCSB CCD API
     (results persisted in db/ccd_smiles_cache.json — reused across runs)
  4. Detect functional groups via FG_SMARTS (substructure match)
  5. Build co-occurrence table:
       rows = amino acid (1-letter code, e.g. H, C, Y)
       cols = functional group name
       values = number of binding events where that FG is present in
                a ligand that contacts that residue

Usage:
    python utils/interaction_analyzer.py
    python utils/interaction_analyzer.py --top 500          # quick test
    python utils/interaction_analyzer.py --local db/BioLiP_nr.txt.gz
    python utils/interaction_analyzer.py --update           # check & re-run if BioLiP changed
"""

import argparse
import gzip
import io
import json
import sys
import time
from collections import defaultdict
from pathlib import Path

import pandas as pd
import requests
from rdkit import Chem, RDLogger

# Suppress RDKit SMILES/SMARTS parse noise — invalid ligands are handled
# explicitly (ligand_to_smiles returns None; detect_fgs returns []).
RDLogger.DisableLog("rdApp.error")
RDLogger.DisableLog("rdApp.warning")

_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT))

from constants.fg_smarts import FG_SMARTS  # noqa: E402

# ── Paths ──────────────────────────────────────────────────────────────────────
BIOLIP_NR_URL    = "https://zhanggroup.org/BioLiP/download/BioLiP_nr.txt.gz"
BIOLIP_NR_CACHE  = _ROOT / "db" / "BioLiP_nr.txt.gz"
OUTPUT_CSV       = _ROOT / "db" / "fg_residue_table.csv"
SMILES_CACHE_PATH = _ROOT / "db" / "ccd_smiles_cache.json"   # persistent SMILES cache
BIOLIP_META_PATH  = _ROOT / "db" / "biolip_metadata.json"    # update-detection metadata

# ── RCSB CCD (Chemical Component Dictionary) API ──────────────────────────────
RCSB_CCD_URL = "https://data.rcsb.org/rest/v1/core/chemcomp/{ccd_id}"
CCD_DELAY    = 0.2   # seconds between RCSB requests

# ── Amino acid 1-letter code → 3-letter name ──────────────────────────────────
AA_1TO3: dict[str, str] = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
    "E": "GLU", "Q": "GLN", "G": "GLY", "H": "HIS", "I": "ILE",
    "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
    "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}
ALL_AA = sorted(AA_1TO3.keys())   # 20 standard amino acids

# ── Ligand exclusion list ──────────────────────────────────────────────────────
# Non-small-molecule ligands that will never have meaningful FG hits
_SKIP_LIGANDS: set[str] = {
    "dna", "rna", "peptide",                    # biopolymers
    "HOH", "DOD",                               # water
    "ZN", "MG", "CA", "FE", "CU", "MN",        # metal ions (2-3 letter)
    "NA", "K", "CL", "BR", "IOD", "PO4", "SO4",
    "GOL", "EDO", "PEG",                        # crystallization additives
}


# ── BioLiP download & update detection ────────────────────────────────────────

def _save_biolip_metadata(dest: Path, content_length: int, etag: str) -> None:
    """Persist download metadata for future update checks."""
    from datetime import date
    meta = {
        "url":            BIOLIP_NR_URL,
        "last_downloaded": str(date.today()),
        "file_size_bytes": dest.stat().st_size,
        "content_length":  content_length,
        "etag":            etag,
    }
    BIOLIP_META_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(BIOLIP_META_PATH, "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2)


def check_biolip_update() -> tuple[bool, str]:
    """Check whether a newer BioLiP_nr.txt.gz is available on the server.

    Uses HTTP HEAD to compare Content-Length and ETag against stored metadata.
    Returns (update_available: bool, description: str).
    """
    if not BIOLIP_META_PATH.exists():
        return True, "No metadata found — treating as first run"

    with open(BIOLIP_META_PATH, encoding="utf-8") as f:
        meta = json.load(f)

    try:
        r = requests.head(BIOLIP_NR_URL, timeout=15)
        server_len  = int(r.headers.get("Content-Length", 0))
        server_etag = r.headers.get("ETag", "").strip('"')
    except Exception as exc:
        return False, f"Cannot reach server: {exc}"

    stored_len  = meta.get("content_length", 0)
    stored_etag = meta.get("etag", "")

    if server_len and server_len != stored_len:
        return True, (
            f"Size changed: stored={stored_len/1e6:.2f} MB  "
            f"server={server_len/1e6:.2f} MB"
        )
    if server_etag and stored_etag and server_etag != stored_etag:
        return True, f"ETag changed: {stored_etag!r} → {server_etag!r}"

    return False, (
        f"Up to date  (size={stored_len/1e6:.2f} MB, "
        f"downloaded={meta.get('last_downloaded', '?')})"
    )


def download_biolip_nr(dest: Path = BIOLIP_NR_CACHE, force: bool = False) -> Path:
    """Download BioLiP_nr.txt.gz in 512 KB chunks and save to dest.

    Skips download if the file already exists, unless force=True.
    Saves metadata to db/biolip_metadata.json after download.
    """
    if dest.exists() and not force:
        print(f"  Using cached: {dest}")
        return dest

    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"  Downloading BioLiP_nr from {BIOLIP_NR_URL} ...")

    r = requests.get(BIOLIP_NR_URL, stream=True, timeout=60)
    r.raise_for_status()

    total = int(r.headers.get("Content-Length", 0))
    downloaded = 0

    with open(dest, "wb") as f:
        for chunk in r.iter_content(chunk_size=524288):   # 512 KB
            f.write(chunk)
            downloaded += len(chunk)
            if total:
                pct = downloaded / total * 100
                print(f"\r  {downloaded/1e6:.1f} / {total/1e6:.1f} MB ({pct:.0f}%)",
                      end="", flush=True)
    print(f"\n  Saved: {dest}")

    _save_biolip_metadata(
        dest,
        content_length=total,
        etag=r.headers.get("ETag", "").strip('"'),
    )
    return dest


# ── Parse ──────────────────────────────────────────────────────────────────────

def _parse_residues(residue_str: str) -> list[str]:
    """Parse BioLiP residue string → list of 1-letter AA codes.

    Input format: "F18 A19 G20 L21 H93" (1-letter + resnum, space-separated)
    Returns only standard amino acid codes (filters out non-standard entries).
    """
    aa_list = []
    for token in residue_str.strip().split():
        if token and token[0].upper() in AA_1TO3:
            aa_list.append(token[0].upper())
    return aa_list


def parse_biolip(filepath: Path, top: int | None = None) -> list[dict]:
    """Parse BioLiP flat file and return list of entry dicts.

    Each dict: {"ligand_id": str, "residues": [str, ...]}

    Args:
        filepath: Path to BioLiP_nr.txt.gz (or uncompressed .txt)
        top:      Limit to first N entries (None = all)
    """
    opener = gzip.open if filepath.suffix == ".gz" else open
    entries: list[dict] = []

    with opener(filepath, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split("\t")
            if len(parts) < 8:
                continue

            ligand_id = parts[4].strip()      # col 5 (0-indexed: 4)
            residue_str = parts[7].strip()    # col 8 (0-indexed: 7)

            if not ligand_id or not residue_str:
                continue
            if ligand_id in _SKIP_LIGANDS:
                continue

            residues = _parse_residues(residue_str)
            if not residues:
                continue

            entries.append({"ligand_id": ligand_id, "residues": residues})

            if top and len(entries) >= top:
                break

    return entries


# ── SMILES lookup — persistent cache ──────────────────────────────────────────

def _load_smiles_cache() -> dict[str, str | None]:
    """Load db/ccd_smiles_cache.json on startup (empty dict if missing)."""
    if SMILES_CACHE_PATH.exists():
        with open(SMILES_CACHE_PATH, encoding="utf-8") as f:
            return json.load(f)
    return {}


def save_smiles_cache() -> None:
    """Persist the in-memory SMILES cache to db/ccd_smiles_cache.json.

    Call this at the end of a run (or periodically) to avoid re-querying
    the same ligands on future runs or after interruptions.
    """
    SMILES_CACHE_PATH.parent.mkdir(parents=True, exist_ok=True)
    tmp = SMILES_CACHE_PATH.with_suffix(".tmp")
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(_smiles_cache, f, indent=2, ensure_ascii=False)
    tmp.replace(SMILES_CACHE_PATH)
    print(f"  SMILES cache saved: {len(_smiles_cache)} entries → {SMILES_CACHE_PATH}")


# Load from disk at import time so every run benefits from previous queries
_smiles_cache: dict[str, str | None] = _load_smiles_cache()


def ligand_to_smiles(ligand_id: str) -> str | None:
    """Fetch SMILES for a PDB 3-letter ligand code via RCSB CCD API.

    Results are persisted in db/ccd_smiles_cache.json across runs.
    Returns None if the ligand is not found or has no valid SMILES.
    """
    if ligand_id in _smiles_cache:
        return _smiles_cache[ligand_id]

    url = RCSB_CCD_URL.format(ccd_id=ligand_id.upper())
    smiles = None
    try:
        r = requests.get(url, timeout=10)
        time.sleep(CCD_DELAY)
        if r.status_code == 200:
            data = r.json()
            # pdbx_chem_comp_descriptor: list of {type, descriptor}
            for desc in data.get("pdbx_chem_comp_descriptor", []):
                if desc.get("type") == "SMILES_CANONICAL":
                    smiles = desc.get("descriptor")
                    break
            # fallback: any SMILES entry
            if not smiles:
                for desc in data.get("pdbx_chem_comp_descriptor", []):
                    if desc.get("type") == "SMILES":
                        smiles = desc.get("descriptor")
                        break
    except Exception:
        smiles = None

    # Validate with RDKit
    if smiles and Chem.MolFromSmiles(smiles) is None:
        smiles = None

    _smiles_cache[ligand_id] = smiles
    return smiles


# ── FG detection ──────────────────────────────────────────────────────────────

_fg_patterns: dict[str, object] = {
    name: Chem.MolFromSmarts(smarts)
    for name, smarts in FG_SMARTS.items()
    if Chem.MolFromSmarts(smarts) is not None
}


def detect_fgs(smiles: str) -> list[str]:
    """Return list of FG names (from FG_SMARTS) present in the given SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []
    return [
        name for name, pattern in _fg_patterns.items()
        if mol.GetSubstructMatches(pattern)
    ]


# ── Table builder ──────────────────────────────────────────────────────────────

def build_interaction_table(entries: list[dict]) -> pd.DataFrame:
    """Build the FG × residue co-occurrence table.

    For each entry: if a ligand contains FG X and contacts residue Y,
    increment table[Y][X] by 1.

    Returns:
        DataFrame with rows = amino acid (1-letter), cols = FG name.
    """
    counts: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    fg_names = list(FG_SMARTS.keys())

    total = len(entries)
    resolved = 0
    skipped_no_smiles = 0
    skipped_no_fg = 0

    unique_ligands = {e["ligand_id"] for e in entries}
    print(f"  {total} binding events, {len(unique_ligands)} unique ligands")
    print("  Resolving SMILES and detecting FGs ...")

    # Pre-fetch SMILES for all unique ligands
    ligand_fgs: dict[str, list[str]] = {}
    for i, lig_id in enumerate(sorted(unique_ligands), 1):
        smiles = ligand_to_smiles(lig_id)
        if smiles is None:
            skipped_no_smiles += 1
            ligand_fgs[lig_id] = []
            continue
        fgs = detect_fgs(smiles)
        ligand_fgs[lig_id] = fgs
        if fgs:
            resolved += 1
        else:
            skipped_no_fg += 1

        if i % 50 == 0:
            print(f"    {i}/{len(unique_ligands)} ligands processed ...", flush=True)

    print(f"  Ligands: {resolved} with FGs, "
          f"{skipped_no_smiles} no SMILES, {skipped_no_fg} no FG match")

    # Build table
    for entry in entries:
        fgs = ligand_fgs.get(entry["ligand_id"], [])
        if not fgs:
            continue
        for aa in set(entry["residues"]):   # deduplicate residues per entry
            for fg in fgs:
                counts[aa][fg] += 1

    # Build DataFrame with all 20 standard AAs as rows
    df = pd.DataFrame(counts).T.fillna(0).astype(int)
    df = df.reindex(index=ALL_AA, columns=fg_names, fill_value=0)
    df.index.name = "residue"

    # Add 3-letter label column
    df.insert(0, "residue_name", df.index.map(AA_1TO3))

    return df


# ── Entry point ────────────────────────────────────────────────────────────────

def main() -> None:
    """Run the full interaction table pipeline."""
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

    parser = argparse.ArgumentParser(description="Build FG × residue interaction table")
    parser.add_argument("--local", default=None,
                        help="Path to local BioLiP_nr.txt.gz (skips download)")
    parser.add_argument("--top", type=int, default=None,
                        help="Limit to first N binding events (for quick testing)")
    parser.add_argument("--output", default=str(OUTPUT_CSV),
                        help="Output CSV path")
    parser.add_argument("--update", action="store_true",
                        help=(
                            "Check if BioLiP_nr has been updated on the server. "
                            "If yes, re-download and rebuild the table. "
                            "SMILES cache is always reused to save API calls."
                        ))
    args = parser.parse_args()

    # ── Update check ──────────────────────────────────────────────────────────
    if args.update:
        print("Checking for BioLiP_nr updates ...")
        has_update, reason = check_biolip_update()
        print(f"  {reason}")
        if not has_update:
            print("Nothing to do.")
            return
        print("  Re-downloading BioLiP_nr ...")
        biolip_path = BIOLIP_NR_CACHE
        download_biolip_nr(biolip_path, force=True)
    else:
        biolip_path = Path(args.local) if args.local else BIOLIP_NR_CACHE
        if not biolip_path.exists():
            print("Downloading BioLiP_nr ...")
            download_biolip_nr(biolip_path)
        else:
            print(f"Using: {biolip_path}")

    # ── Parse ─────────────────────────────────────────────────────────────────
    print(f"Parsing BioLiP (top={args.top or 'all'}) ...")
    entries = parse_biolip(biolip_path, top=args.top)
    print(f"  Parsed {len(entries)} binding events")
    print(f"  SMILES cache pre-loaded: {len(_smiles_cache)} entries")

    # ── Build table ───────────────────────────────────────────────────────────
    print("Building interaction table ...")
    df = build_interaction_table(entries)

    # ── Save table ────────────────────────────────────────────────────────────
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=True)
    print(f"\nSaved: {out}")
    print()
    print(df.to_string())

    # ── Persist SMILES cache for future runs ──────────────────────────────────
    save_smiles_cache()


if __name__ == "__main__":
    main()
