"""PDB batch fetcher: search protein family → download → extract FG×residue interactions.

Builds/updates  db/fg_residue_map.json  organised by protein family name.
Supports resume from breakpoint (tracks processed PDB IDs between runs).

Requirements
------------
  biopython  (conda install -c conda-forge biopython)

Usage
-----
  # First run — search + process first 100 entries
  python utils/pdb_batch_fetcher.py --family "serine protease" --limit 100

  # Resume — skips already-done entries automatically
  python utils/pdb_batch_fetcher.py --family "serine protease" --limit 100

  # Pfam-precise search (instead of full-text)
  python utils/pdb_batch_fetcher.py --family kinase --pfam PF00069

  # Dry-run: search only, no download
  python utils/pdb_batch_fetcher.py --family "nuclear receptor" --dry-run

db/fg_residue_map.json schema
-------------------------------
{
  "last_updated": "YYYY-MM-DD",
  "families": {
    "<family_key>": {
      "query":           "<original query string>",
      "total_found":     <int>,
      "pdb_ids_all":     ["1ABC", ...],    // full search result (resume source)
      "pdb_ids_done":    ["1ABC", ...],    // all processed entries (won't be retried)
      "pdb_ids_failed":  {"2DEF": "download_failed", ...}, // hard failures (retryable)
      "fg_residue_counts": {
        "<FG name>": {"<AA 1-letter>": <int>, ...},
        ...
      }
    }
  }
}

Hard failures (download_failed, parse_error) are NOT added to pdb_ids_done,
so they will be retried on the next run.
Soft skips (no_ligands, no_fg_match) ARE added to pdb_ids_done with a note
in pdb_ids_failed — they won't be re-downloaded.
"""

import argparse
import json
import sys
import time
from datetime import date
from pathlib import Path
from typing import Optional

import requests
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.error")
RDLogger.DisableLog("rdApp.warning")

_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT))

from constants.fg_smarts import FG_SMARTS  # noqa: E402

# ── Paths ──────────────────────────────────────────────────────────────────────
TMP_DIR  = _ROOT / "tmp"
MAP_PATH = _ROOT / "db" / "fg_residue_map.json"

# ── RCSB API endpoints ─────────────────────────────────────────────────────────
RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_PDB_URL    = "https://files.rcsb.org/download/{pdb_id}.pdb"
RCSB_CCD_URL    = "https://data.rcsb.org/rest/v1/core/chemcomp/{ccd_id}"

API_DELAY      = 0.15   # seconds between RCSB requests (be polite)
CONTACT_CUTOFF = 4.5    # Angstroms — standard protein–ligand contact distance
SEARCH_ROWS    = 200    # results per RCSB search page
CHECKPOINT     = 25     # save progress JSON every N entries

# ── Residue / ligand exclusion lists ──────────────────────────────────────────
SKIP_LIGANDS: frozenset[str] = frozenset({
    # Water & crystallisation solvents
    "HOH", "DOD", "WAT", "H2O",
    "GOL", "EDO", "PEG", "IPA", "DMS", "MPD", "HED", "BME",
    # Buffers / common additives
    "SO4", "PO4", "ACT", "ACE", "FMT", "TRS", "EPE", "MES", "CIT",
    # Metal ions (2–3 letter CCD codes)
    "NA", "K", "ZN", "MG", "CA", "FE", "CU", "MN", "CO", "NI",
    "CD", "HG", "PT", "AU", "AG", "LI", "SR", "BA",
    # Halogens / simple anions
    "CL", "BR", "IOD", "F", "FLU",
    # Nucleotides (DNA / RNA)
    "DA", "DT", "DG", "DC", "DI", "DU",
    "AMP", "GMP", "CMP", "UMP", "ATP", "GTP", "ADP", "GDP",
    # Modified amino acids / linkers
    "MSE", "SEC", "PYL",
})

# ── Amino acid mappings ────────────────────────────────────────────────────────
_AA_1TO3: dict[str, str] = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
    "E": "GLU", "Q": "GLN", "G": "GLY", "H": "HIS", "I": "ILE",
    "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
    "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}
AA_3TO1: dict[str, str] = {v: k for k, v in _AA_1TO3.items()}

# ── Pre-compile SMARTS once at import time ─────────────────────────────────────
_SMARTS_PATTERNS: dict[str, object] = {
    name: Chem.MolFromSmarts(smarts)
    for name, smarts in FG_SMARTS.items()
    if Chem.MolFromSmarts(smarts) is not None
}

# ── In-memory SMILES cache (session scope) ─────────────────────────────────────
_smiles_cache: dict[str, Optional[str]] = {}


# ═══════════════════════════════════════════════════════════════════════════════
# 1. RCSB Search
# ═══════════════════════════════════════════════════════════════════════════════

def search_pdb_family(
    family: str,
    pfam_id: Optional[str] = None,
    max_results: int = 5000,
) -> list[str]:
    """Search RCSB for PDB entry IDs that match a protein family.

    Args:
        family:      Family name used in full-text search (all RCSB text fields).
        pfam_id:     Pfam accession (e.g. "PF00069"); when given, runs an exact
                     annotation-match query instead of full-text search.
        max_results: Hard cap on how many PDB IDs to return.

    Returns:
        List of PDB IDs (uppercase 4-char codes).
    """
    if pfam_id:
        query_node = {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_entity_annotation.annotation_lineage.id",
                "operator": "exact_match",
                "value": pfam_id.upper(),
            },
        }
        label = f"Pfam:{pfam_id.upper()}"
    else:
        query_node = {
            "type": "terminal",
            "service": "full_text",
            "parameters": {"value": family},
        }
        label = f'"{family}"'

    pdb_ids: list[str] = []
    start = 0

    print(f"  Searching RCSB for {label} ...", flush=True)

    while len(pdb_ids) < max_results:
        payload = {
            "query": query_node,
            "return_type": "entry",
            "request_options": {
                "paginate": {"start": start, "rows": SEARCH_ROWS}
            },
        }
        try:
            r = requests.post(RCSB_SEARCH_URL, json=payload, timeout=30)
            r.raise_for_status()
            data = r.json()
        except Exception as exc:
            print(f"\n  ! Search error at offset {start}: {exc}", file=sys.stderr)
            break

        total     = data.get("total_count", 0)
        result_set = data.get("result_set", [])
        if not result_set:
            break

        for entry in result_set:
            pdb_ids.append(entry["identifier"].upper())

        cap = min(total, max_results)
        print(f"    {len(pdb_ids):>5} / {cap} retrieved ...", end="\r", flush=True)

        start += SEARCH_ROWS
        if start >= total or len(pdb_ids) >= max_results:
            break
        time.sleep(API_DELAY)

    print(f"\n  Total found: {len(pdb_ids)} entries for {label}")
    return pdb_ids[:max_results]


# ═══════════════════════════════════════════════════════════════════════════════
# 2. PDB download
# ═══════════════════════════════════════════════════════════════════════════════

def download_pdb(pdb_id: str, dest_dir: Path) -> Optional[Path]:
    """Download the PDB-format file for pdb_id into dest_dir.

    Skips download if the file already exists (safe restart).
    Returns the path to the saved file, or None on failure.
    """
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / f"{pdb_id.lower()}.pdb"
    if dest.exists():
        return dest

    url = RCSB_PDB_URL.format(pdb_id=pdb_id.lower())
    for attempt in range(1, 3):
        try:
            r = requests.get(url, timeout=30, stream=True)
            if r.status_code == 404:
                return None   # structure not available in PDB format
            r.raise_for_status()
            with open(dest, "wb") as f:
                for chunk in r.iter_content(chunk_size=131072):  # 128 KB
                    f.write(chunk)
            time.sleep(API_DELAY)
            return dest
        except Exception as exc:
            if attempt == 2:
                print(f"    ! Download failed ({pdb_id}): {exc}", file=sys.stderr)
                return None
            time.sleep(1.5)
    return None


# ═══════════════════════════════════════════════════════════════════════════════
# 3. Chemistry helpers
# ═══════════════════════════════════════════════════════════════════════════════

def ligand_to_smiles(ccd_id: str) -> Optional[str]:
    """Fetch canonical SMILES for a PDB CCD ligand code via RCSB REST API.

    Results are cached in memory for the duration of the process.
    Returns None if the ligand is unknown or has no valid SMILES.
    """
    if ccd_id in _smiles_cache:
        return _smiles_cache[ccd_id]

    url = RCSB_CCD_URL.format(ccd_id=ccd_id.upper())
    smiles: Optional[str] = None
    try:
        r = requests.get(url, timeout=10)
        time.sleep(API_DELAY)
        if r.status_code == 200:
            data = r.json()
            descriptors = data.get("pdbx_chem_comp_descriptor", [])
            # Prefer SMILES_CANONICAL, fall back to any SMILES type
            for desc in descriptors:
                if desc.get("type") == "SMILES_CANONICAL":
                    smiles = desc.get("descriptor")
                    break
            if not smiles:
                for desc in descriptors:
                    if "SMILES" in desc.get("type", ""):
                        smiles = desc.get("descriptor")
                        break
    except Exception:
        pass

    # Validate with RDKit (silently discards unparsable SMILES)
    if smiles and Chem.MolFromSmiles(smiles) is None:
        smiles = None

    _smiles_cache[ccd_id] = smiles
    return smiles


def detect_fgs(smiles: str) -> list[str]:
    """Return FG names (from FG_SMARTS) found via substructure match in smiles."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []
    return [
        name for name, pat in _SMARTS_PATTERNS.items()
        if pat is not None and mol.GetSubstructMatches(pat)
    ]


# ═══════════════════════════════════════════════════════════════════════════════
# 4. Structure analysis (BioPython)
# ═══════════════════════════════════════════════════════════════════════════════

def _require_biopython() -> None:
    """Raise ImportError with install hint if BioPython is not available."""
    try:
        import Bio  # noqa: F401
    except ImportError:
        raise ImportError(
            "BioPython is required for PDB structure analysis.\n"
            "Install with:\n"
            "  conda install -c conda-forge biopython\n"
            "  # or: pip install biopython"
        ) from None


def extract_contacts(pdb_path: Path) -> list[tuple[str, list[str]]]:
    """Find small-molecule ligand ↔ protein-residue contacts in a PDB file.

    Algorithm:
        1. Build a NeighborSearch index over all protein ATOM-record atoms.
        2. For each HETATM ligand (not in SKIP_LIGANDS), query the index
           for all protein residues within CONTACT_CUTOFF Angstroms of any
           ligand atom.
        3. Return the 3-letter CCD code and the set of contacting AA (1-letter).

    Args:
        pdb_path: Path to the PDB-format file.

    Returns:
        List of (ccd_id, [aa_1letter, ...]) — one tuple per unique ligand instance.
        Raises ValueError if BioPython cannot parse the file.
    """
    from Bio.PDB import PDBParser
    from Bio.PDB.NeighborSearch import NeighborSearch

    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("s", str(pdb_path))
    except Exception as exc:
        raise ValueError(f"PDBParser failed on {pdb_path.name}: {exc}") from exc

    # Index all protein (ATOM-record) atoms for O(log N) spatial lookup
    protein_atoms = [
        atom for atom in structure.get_atoms()
        if atom.get_parent().id[0] == " "   # standard residue flag
    ]
    if not protein_atoms:
        return []

    ns = NeighborSearch(protein_atoms)
    results: list[tuple[str, list[str]]] = []
    seen_ligands: set[str] = set()

    for model in structure:
        for chain in model:
            for residue in chain:
                het_flag = residue.id[0]

                # BioPython het_flag convention:
                #   " "      → standard residue (ATOM record)
                #   "W"      → water
                #   "H_XXX"  → HETATM with residue name XXX
                if not het_flag.startswith("H_"):
                    continue

                ccd_id = residue.get_resname().strip().upper()
                if ccd_id in SKIP_LIGANDS:
                    continue

                # Deduplicate: same ligand may appear in multiple models
                lig_key = f"{chain.id}:{residue.id[1]}:{ccd_id}"
                if lig_key in seen_ligands:
                    continue
                seen_ligands.add(lig_key)

                # Spatial query: find protein residues within cutoff
                nearby: set[str] = set()
                for atom in residue:
                    for nb_res in ns.search(
                        atom.get_coord(), CONTACT_CUTOFF, level="R"
                    ):
                        aa1 = AA_3TO1.get(nb_res.get_resname().strip())
                        if aa1:
                            nearby.add(aa1)

                if nearby:
                    results.append((ccd_id, sorted(nearby)))

    return results


# ═══════════════════════════════════════════════════════════════════════════════
# 5. Map JSON I/O
# ═══════════════════════════════════════════════════════════════════════════════

def _empty_family(query: str) -> dict:
    """Return a freshly initialised family entry."""
    return {
        "query":            query,
        "total_found":      0,
        "pdb_ids_all":      [],
        "pdb_ids_done":     [],
        "pdb_ids_failed":   {},   # {pdb_id: reason_string}
        "fg_residue_counts": {},  # {fg_name: {aa_1letter: int}}
    }


def load_map(path: Path = MAP_PATH) -> dict:
    """Load db/fg_residue_map.json, or return an empty skeleton."""
    if path.exists():
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    return {"last_updated": str(date.today()), "families": {}}


def save_map(data: dict, path: Path = MAP_PATH) -> None:
    """Atomically write fg_residue_map.json (write *.tmp → rename)."""
    data["last_updated"] = str(date.today())
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    tmp.replace(path)  # atomic on POSIX; near-atomic on Windows


# ═══════════════════════════════════════════════════════════════════════════════
# 6. Core processor (one PDB entry)
# ═══════════════════════════════════════════════════════════════════════════════

def process_entry(pdb_id: str, fam: dict) -> str:
    """Download and analyse one PDB entry; update fam['fg_residue_counts'].

    The temporary PDB file is always deleted after analysis (success or failure).

    Returns one of:
        "ok"               interactions found and counts updated
        "no_ligands"       structure has no eligible HETATM ligands
        "no_fg_match"      ligands found but none matched FG_SMARTS
        "download_failed"  RCSB file unavailable or network error
        "parse_error:..."  BioPython could not parse the file
    """
    # ── Download ───────────────────────────────────────────────────────────────
    pdb_path = download_pdb(pdb_id, TMP_DIR)
    if pdb_path is None:
        return "download_failed"

    # ── Analyse (always clean up tmp file) ────────────────────────────────────
    try:
        contacts = extract_contacts(pdb_path)
    except ValueError as exc:
        return f"parse_error: {exc}"
    finally:
        if pdb_path.exists():
            pdb_path.unlink()

    if not contacts:
        return "no_ligands"

    # ── Accumulate FG × residue counts ────────────────────────────────────────
    fg_counts: dict = fam["fg_residue_counts"]
    any_hit = False

    for ccd_id, residues in contacts:
        smiles = ligand_to_smiles(ccd_id)
        if not smiles:
            continue
        fgs = detect_fgs(smiles)
        if not fgs:
            continue

        any_hit = True
        for fg in fgs:
            aa_map = fg_counts.setdefault(fg, {})
            for aa in residues:
                aa_map[aa] = aa_map.get(aa, 0) + 1

    return "ok" if any_hit else "no_fg_match"


# ═══════════════════════════════════════════════════════════════════════════════
# 7. Entry point
# ═══════════════════════════════════════════════════════════════════════════════

def main() -> None:
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description=(
            "Search PDB by protein family, download structures, "
            "and build a FG×residue interaction map (db/fg_residue_map.json)."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--family", required=True,
        help='Protein family name, e.g. "serine protease", "kinase"',
    )
    parser.add_argument(
        "--pfam", default=None,
        help="Pfam accession for precise search, e.g. PF00069 (optional)",
    )
    parser.add_argument(
        "--limit", type=int, default=None,
        help="Max entries to process this run (incremental; omit to process all)",
    )
    parser.add_argument(
        "--max-results", type=int, default=5000,
        help="Max PDB IDs to retrieve from RCSB search",
    )
    parser.add_argument(
        "--cutoff", type=float, default=CONTACT_CUTOFF,
        help="Ligand–residue contact distance cutoff in Angstroms",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Search only — report counts without downloading anything",
    )
    parser.add_argument(
        "--output", default=str(MAP_PATH),
        help="Output JSON path",
    )
    args = parser.parse_args()

    # Allow CLI override of global cutoff
    global CONTACT_CUTOFF
    CONTACT_CUTOFF = args.cutoff

    map_path   = Path(args.output)
    family_key = args.family.lower().strip()

    # ── Load / initialise the map ──────────────────────────────────────────────
    data = load_map(map_path)
    if family_key not in data["families"]:
        data["families"][family_key] = _empty_family(args.family)
    fam = data["families"][family_key]

    # ── Step 1: Search (skip if ID list already cached) ───────────────────────
    if not fam["pdb_ids_all"]:
        pdb_ids = search_pdb_family(
            args.family,
            pfam_id=args.pfam,
            max_results=args.max_results,
        )
        fam["pdb_ids_all"] = pdb_ids
        fam["total_found"] = len(pdb_ids)
        save_map(data, map_path)
    else:
        n_all  = len(fam["pdb_ids_all"])
        n_done = len(fam["pdb_ids_done"])
        n_fail = len(fam["pdb_ids_failed"])
        print(
            f"  Resuming [{family_key}]:  "
            f"total={n_all}  done={n_done}  failed={n_fail}  "
            f"remaining={n_all - n_done - n_fail}"
        )

    if args.dry_run:
        n_all  = len(fam["pdb_ids_all"])
        n_done = len(fam["pdb_ids_done"])
        n_fail = len(fam["pdb_ids_failed"])
        print(
            f"  [dry-run]  total={n_all}  done={n_done}  "
            f"failed={n_fail}  remaining={n_all - n_done - n_fail}"
        )
        return

    # ── Step 2: Determine remaining entries ───────────────────────────────────
    _require_biopython()

    done_set   = set(fam["pdb_ids_done"])
    failed_set = set(fam["pdb_ids_failed"].keys())
    remaining  = [
        pid for pid in fam["pdb_ids_all"]
        if pid not in done_set and pid not in failed_set
    ]
    if args.limit:
        remaining = remaining[: args.limit]

    n_remaining = len(remaining)
    if not n_remaining:
        print("  Nothing left to process. All entries are done!")
        return

    print(
        f"  Processing {n_remaining} entries  "
        f"(done={len(done_set)}, hard-failed={len(failed_set)}) ..."
    )
    TMP_DIR.mkdir(parents=True, exist_ok=True)

    # ── Step 3: Process loop ──────────────────────────────────────────────────
    ok_count   = 0
    skip_count = 0
    fail_count = 0

    for i, pdb_id in enumerate(remaining, 1):
        status = process_entry(pdb_id, fam)

        if status == "ok":
            fam["pdb_ids_done"].append(pdb_id)
            ok_count += 1
            marker = "✓"

        elif status in ("no_ligands", "no_fg_match"):
            # Processed but nothing useful → mark done so we don't retry
            fam["pdb_ids_done"].append(pdb_id)
            fam["pdb_ids_failed"][pdb_id] = status   # informational only
            skip_count += 1
            marker = "–"

        else:
            # Hard failure (download_failed / parse_error) → NOT done → retryable
            fam["pdb_ids_failed"][pdb_id] = status
            fail_count += 1
            marker = "✗"

        print(
            f"  [{i:>5}/{n_remaining}] {pdb_id}  {status:<35} {marker}",
            flush=True,
        )

        # Periodic checkpoint (atomic save)
        if i % CHECKPOINT == 0:
            save_map(data, map_path)
            print(
                f"  ── checkpoint @ {i} ──  "
                f"ok={ok_count}  skip={skip_count}  fail={fail_count}",
                flush=True,
            )

    # ── Final save ────────────────────────────────────────────────────────────
    save_map(data, map_path)

    # Remove tmp dir if now empty
    try:
        TMP_DIR.rmdir()
    except OSError:
        pass

    # ── Summary ───────────────────────────────────────────────────────────────
    width = 52
    print(f"\n{'='*width}")
    print(f"  Family     : {family_key}")
    print(f"  Processed  : {ok_count + skip_count + fail_count}")
    print(f"    ✓  ok (FGs found)  : {ok_count:>5}")
    print(f"    –  no useful data  : {skip_count:>5}")
    print(f"    ✗  hard failures   : {fail_count:>5}")
    print(f"  FG types recorded  : {len(fam['fg_residue_counts'])}")
    print(f"  Output     : {map_path}")
    print(f"{'='*width}")


if __name__ == "__main__":
    main()
