"""Step 2b — Extract 3D residue-ligand interaction poses from BioLiP structures.

For each FG × residue pair (top-ranked by fg_residue_table.csv):
  1. Re-parse BioLiP with full metadata (PDB ID, chain, resolution, lig serial)
  2. Select representative structures (filtered by resolution, capped per pair)
  3. Download PDB structures from RCSB
  4. Extract 3D data with BioPython
  5. Save three complementary outputs:

     db/residue_3d_poses.json
         Cα coordinate + ligand centroid + distance per binding event.
         Primary input for Step 3 (structural motif search).

     db/local_env/{fg_safe}_{residue}.sdf
         Multi-molecule SDF: one mol = one ligand (all HETATM atoms).
         Useful for visualisation, manual inspection, and MCS-based alignment.

     db/pharmacophore_stats.json
         Per (FG, residue) distance/frequency statistics.
         Used as a confidence filter in the final prediction report.

Usage
-----
    python utils/pose_extractor.py
    python utils/pose_extractor.py --top-pairs 15 --max-per-pair 25
    python utils/pose_extractor.py --resolution 2.0
    python utils/pose_extractor.py --fg "Phenol" --residue Y
    python utils/pose_extractor.py --dry-run   # show plan, no downloads
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import re
import shutil
import sys
import time
from collections import defaultdict
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterator

import numpy as np
import pandas as pd
import requests
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.error")
RDLogger.DisableLog("rdApp.warning")

_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT))

from constants.fg_smarts import FG_SMARTS  # noqa: E402

# ── Paths ──────────────────────────────────────────────────────────────────────
BIOLIP_NR_CACHE   = _ROOT / "db" / "BioLiP_nr.txt.gz"
SMILES_CACHE_PATH = _ROOT / "db" / "ccd_smiles_cache.json"
FG_RESIDUE_TABLE  = _ROOT / "db" / "fg_residue_table.csv"

POSES_JSON        = _ROOT / "db" / "residue_3d_poses.json"
STATS_JSON        = _ROOT / "db" / "pharmacophore_stats.json"
LOCAL_ENV_DIR     = _ROOT / "db" / "local_env"

TMP_DIR           = _ROOT / "tmp"
RCSB_PDB_URL      = "https://files.rcsb.org/download/{pdb_id}.pdb"

# ── Tunable constants ─────────────────────────────────────────────────────────
DEFAULT_MAX_PER_PAIR = 30     # max structures downloaded per FG × residue pair
DEFAULT_TOP_PAIRS    = 20     # number of top FG × residue pairs to process
DEFAULT_MAX_RESOL    = 2.5    # max resolution (Å) to accept; NMR (None) accepted
CONTACT_CUTOFF       = 4.5    # Å — consistent with pdb_batch_fetcher
PDB_DELAY            = 0.15   # seconds between RCSB downloads

# ── BioLiP column indices ─────────────────────────────────────────────────────
_COL_PDBID   = 0
_COL_CHAIN   = 1   # receptor (protein) chain
_COL_RESOL   = 2   # resolution in Å (empty string for NMR)
_COL_SITE    = 3   # binding site ID within PDB
_COL_LIGID   = 4   # ligand CCD code
_COL_LGCHN   = 5   # ligand chain
_COL_LGSER   = 6   # ligand residue sequence number (author)
_COL_RESIDS  = 7   # binding residues "AA1+resnum ..." author numbering
_COL_ECNUM   = 11  # EC number (may be empty)
_COL_UNIPROT = 17  # UniProt accession (may be empty)

# Non-small-molecule ligands — same list as interaction_analyzer.py
_SKIP_LIGANDS: set[str] = {
    "dna", "rna", "peptide",
    "HOH", "DOD",
    "ZN", "MG", "CA", "FE", "CU", "MN",
    "NA", "K", "CL", "BR", "IOD", "PO4", "SO4",
    "GOL", "EDO", "PEG",
}

# Standard amino acid 1-letter codes
_AA_CODES: set[str] = set("ARNDCQEGHILKMFPSTWYV")

# 3-letter → 1-letter reverse mapping (for BioPython residue name lookup)
_AA_3TO1: dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


# ── Data containers ───────────────────────────────────────────────────────────

@dataclass
class BioLiPEntry:
    """One binding event from BioLiP, with full metadata for 3D extraction."""
    pdb_id: str
    rec_chain: str
    resolution: float | None   # None for NMR
    site_id: str
    lig_id: str
    lig_chain: str
    lig_serial: int
    res_tokens: list[str]       # ["F18", "A19", ...] author-numbered
    ec_number: str
    uniprot_id: str


@dataclass
class PoseRecord:
    """3D contact data extracted from one PDB binding event."""
    pdb_id: str
    lig_id: str
    lig_chain: str
    lig_serial: int             # ligand residue seq number in PDB (author)
    rec_chain: str
    res_code: str               # 1-letter AA
    res_num: int                # protein residue seq number
    fg_name: str
    ca_coord: list[float] | None
    lig_centroid: list[float]
    ca_lig_dist: float | None   # Cα to ligand centroid
    min_atom_dist: float        # min any-atom to any-atom distance
    resolution: float | None
    ec_number: str


# ── BioLiP parsing (full columns) ─────────────────────────────────────────────

def _parse_resolution(s: str) -> float | None:
    """Convert resolution string to float; return None for NMR / missing."""
    s = s.strip()
    if not s or s.lower() in ("", "nmr", "na", "null", "0.0", "0"):
        return None
    try:
        v = float(s)
        return v if v > 0 else None
    except ValueError:
        return None


def _parse_res_tokens(residue_str: str) -> list[str]:
    """Parse 'F18 A19 G20 ...' → ['F18', 'A19', 'G20', ...]

    Keeps only tokens whose first character is a standard AA 1-letter code.
    """
    tokens = []
    for tok in residue_str.strip().split():
        if tok and tok[0].upper() in _AA_CODES:
            tokens.append(tok)
    return tokens


def parse_biolip_full(filepath: Path, top: int | None = None) -> list[BioLiPEntry]:
    """Parse BioLiP_nr.txt.gz and return full entry list.

    Unlike interaction_analyzer.parse_biolip, this version keeps all metadata
    needed for 3D structure extraction (PDB ID, chains, serial numbers,
    resolution, UniProt).

    Args:
        filepath: Path to BioLiP_nr.txt.gz (or .txt)
        top:      Parse at most this many entries (None = all)
    """
    opener = gzip.open if str(filepath).endswith(".gz") else open
    entries: list[BioLiPEntry] = []

    with opener(filepath, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 8:
                continue

            lig_id = parts[_COL_LIGID].strip()
            if not lig_id or lig_id in _SKIP_LIGANDS:
                continue

            res_str = parts[_COL_RESIDS].strip()
            res_tokens = _parse_res_tokens(res_str)
            if not res_tokens:
                continue

            try:
                lig_serial = int(parts[_COL_LGSER].strip())
            except (ValueError, IndexError):
                lig_serial = 0

            entry = BioLiPEntry(
                pdb_id    = parts[_COL_PDBID].strip().lower(),
                rec_chain = parts[_COL_CHAIN].strip(),
                resolution= _parse_resolution(
                    parts[_COL_RESOL] if len(parts) > _COL_RESOL else ""),
                site_id   = parts[_COL_SITE].strip() if len(parts) > _COL_SITE else "",
                lig_id    = lig_id,
                lig_chain = parts[_COL_LGCHN].strip() if len(parts) > _COL_LGCHN else "",
                lig_serial= lig_serial,
                res_tokens= res_tokens,
                ec_number = parts[_COL_ECNUM].strip() if len(parts) > _COL_ECNUM else "",
                uniprot_id= parts[_COL_UNIPROT].strip() if len(parts) > _COL_UNIPROT else "",
            )
            entries.append(entry)

            if top and len(entries) >= top:
                break

    return entries


# ── FG detection (shared with interaction_analyzer) ───────────────────────────

_SMARTS_PATTERNS: dict[str, Chem.Mol] = {
    name: Chem.MolFromSmarts(s)
    for name, s in FG_SMARTS.items()
    if Chem.MolFromSmarts(s) is not None
}


def _detect_fgs(smiles: str) -> list[str]:
    """Return FG names present in the molecule (subset of FG_SMARTS keys)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []
    return [
        name for name, pat in _SMARTS_PATTERNS.items()
        if pat is not None and mol.GetSubstructMatches(pat)
    ]


def _load_smiles_cache() -> dict[str, str | None]:
    """Load ccd_smiles_cache.json (built by interaction_analyzer.py)."""
    if SMILES_CACHE_PATH.exists():
        with open(SMILES_CACHE_PATH, encoding="utf-8") as f:
            return json.load(f)
    return {}


# ── Index: FG × residue → [BioLiPEntry] ──────────────────────────────────────

def build_fg_residue_index(
    entries: list[BioLiPEntry],
    smiles_cache: dict[str, str | None],
) -> dict[str, dict[str, list[BioLiPEntry]]]:
    """Build {fg_name: {res_1letter: [BioLiPEntry, ...]}} lookup.

    Iterates entries, checks if the ligand SMILES is known, detects FGs,
    and records each entry under all FG × residue combinations it contributes.
    """
    index: dict[str, dict[str, list[BioLiPEntry]]] = defaultdict(
        lambda: defaultdict(list)
    )

    n_no_smiles = n_no_fg = n_ok = 0

    for entry in entries:
        smiles = smiles_cache.get(entry.lig_id)
        if not smiles:
            n_no_smiles += 1
            continue

        fgs = _detect_fgs(smiles)
        if not fgs:
            n_no_fg += 1
            continue

        for fg in fgs:
            unique_res = {tok[0].upper() for tok in entry.res_tokens}
            for res in unique_res:
                index[fg][res].append(entry)

        n_ok += 1

    print(f"  Index built: {n_ok} entries with FGs, "
          f"{n_no_smiles} no SMILES, {n_no_fg} no FG match")
    return index


# ── Representative selection ──────────────────────────────────────────────────

def select_representatives(
    index: dict[str, dict[str, list[BioLiPEntry]]],
    top_pairs: int,
    max_per_pair: int,
    max_resolution: float | None,
    fg_filter: str | None,
    res_filter: str | None,
    top_per_fg: int | None = None,
) -> dict[str, dict[str, list[BioLiPEntry]]]:
    """Select top representative entries per FG × residue pair.

    Two selection modes (mutually exclusive):

    --top-per-fg N  (recommended for comprehensive coverage)
        For each FG, pick the top N residue partners by BioLiP count.
        Guarantees all FG types are represented even if they are rare.
        Example: --top-per-fg 3  →  21 FGs × 3 residues = up to 63 pairs.

    --top-pairs N  (default, global ranking)
        Pick the globally top N FG × residue pairs by BioLiP count.
        May focus entirely on the most common FGs (Hydroxyl, Ether) if N is small.

    Within each pair, structures are sorted by resolution (best first),
    deduplicated by PDB ID, and capped at max_per_pair.
    """
    # Load reference table for ranking
    if FG_RESIDUE_TABLE.exists():
        ref_df = pd.read_csv(FG_RESIDUE_TABLE, index_col=0)
    else:
        ref_df = None

    def _pair_score(fg: str, res: str) -> int:
        """Higher = more common interaction in BioLiP."""
        if ref_df is None:
            return len(index.get(fg, {}).get(res, []))
        try:
            return int(ref_df.at[res, fg])
        except (KeyError, TypeError):
            return 0

    # Enumerate all valid pairs, apply filters
    all_pairs: list[tuple[str, str, int]] = []
    for fg, res_dict in index.items():
        if fg_filter and fg != fg_filter:
            continue
        for res, ents in res_dict.items():
            if res_filter and res != res_filter.upper():
                continue
            all_pairs.append((fg, res, _pair_score(fg, res)))

    # ── Mode 1: top N per FG (ensures coverage of all FG types) ──────────────
    if top_per_fg is not None:
        from itertools import groupby
        # Group pairs by FG, then select top top_per_fg within each FG
        by_fg: dict[str, list[tuple[str, str, int]]] = defaultdict(list)
        for fg, res, score in all_pairs:
            by_fg[fg].append((fg, res, score))

        selected_pairs = []
        for fg, pairs in sorted(by_fg.items()):
            pairs.sort(key=lambda x: x[2], reverse=True)
            selected_pairs.extend(pairs[:top_per_fg])

        print(f"  Mode: top-{top_per_fg}-per-FG "
              f"→ {len(selected_pairs)} pairs across {len(by_fg)} FGs")

    # ── Mode 2: global top N pairs ────────────────────────────────────────────
    else:
        all_pairs.sort(key=lambda x: x[2], reverse=True)
        selected_pairs = all_pairs[:top_pairs]
        print(f"  Mode: top-{top_pairs} global pairs "
              f"(of {len(all_pairs)} total)")

    result: dict[str, dict[str, list[BioLiPEntry]]] = defaultdict(dict)

    for fg, res, score in selected_pairs:
        candidates = index[fg][res]

        # Filter by resolution
        if max_resolution is not None:
            candidates = [
                e for e in candidates
                if e.resolution is None or e.resolution <= max_resolution
            ]

        # Sort: lower resolution = better quality
        def _sort_key(e: BioLiPEntry) -> float:
            return e.resolution if e.resolution is not None else 99.0

        candidates.sort(key=_sort_key)

        # Deduplicate by PDB ID
        seen_pdb: set[str] = set()
        unique: list[BioLiPEntry] = []
        for e in candidates:
            if e.pdb_id not in seen_pdb:
                seen_pdb.add(e.pdb_id)
                unique.append(e)
            if len(unique) >= max_per_pair:
                break

        result[fg][res] = unique
        print(f"    {fg:26s} x {res}: {score:6,d} events -> {len(unique):3d} representatives")

    return result


# ── PDB download ──────────────────────────────────────────────────────────────

def _pdb_cache_path(pdb_id: str) -> Path:
    return TMP_DIR / f"{pdb_id.lower()}.pdb"


def download_pdb(pdb_id: str, delay: float = PDB_DELAY) -> Path | None:
    """Download PDB file to tmp/; return path or None on failure.

    Re-uses cached file if already present.
    """
    TMP_DIR.mkdir(parents=True, exist_ok=True)
    dest = _pdb_cache_path(pdb_id)
    if dest.exists():
        return dest

    url = RCSB_PDB_URL.format(pdb_id=pdb_id.lower())
    try:
        r = requests.get(url, timeout=30)
        if r.status_code == 200:
            dest.write_bytes(r.content)
            time.sleep(delay)
            return dest
        else:
            return None
    except Exception:
        return None


# ── 3D extraction ──────────────────────────────────────────────────────────────

def _vec_dist(a: list[float], b: list[float]) -> float:
    """Euclidean distance between two 3D points."""
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def extract_pose_records(
    pdb_path: Path,
    entry: BioLiPEntry,
    fg_name: str,
    cutoff: float = CONTACT_CUTOFF,
) -> list[PoseRecord]:
    """Extract PoseRecords for one BioLiP entry from the PDB file.

    Steps:
      1. Find the ligand HETATM residue (lig_chain + lig_serial).
      2. Compute ligand centroid from all heavy atoms.
      3. For each residue token in entry.res_tokens:
           a. Locate the residue in the protein chain.
           b. Get Cα coordinate.
           c. Compute distances (Cα→centroid, min atom-to-atom).
      4. Return one PoseRecord per (fg_name, residue) combination.

    Returns empty list if the ligand cannot be found.
    """
    try:
        from Bio.PDB import PDBParser
        from Bio.PDB.Residue import Residue as BioResidue
    except ImportError:
        raise ImportError("BioPython is required. Install: conda install -c conda-forge biopython")

    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("x", str(pdb_path))
    except Exception:
        return []

    model = structure[0]

    # ── Find ligand residue ────────────────────────────────────────────────────
    lig_res = None
    # Primary: match lig_chain + lig_serial + HETATM
    for chain in model:
        if chain.id != entry.lig_chain:
            continue
        for res in chain:
            if res.id[1] == entry.lig_serial and res.id[0].startswith("H_"):
                lig_res = res
                break
        if lig_res:
            break
    # Fallback A: any chain, same lig_serial
    if lig_res is None:
        for chain in model:
            for res in chain:
                if res.id[1] == entry.lig_serial and res.id[0].startswith("H_"):
                    lig_res = res
                    break
            if lig_res:
                break
    # Fallback B: first HETATM matching lig_id by name
    if lig_res is None:
        for chain in model:
            for res in chain:
                if res.id[0].startswith("H_") and res.resname.strip() == entry.lig_id.strip():
                    lig_res = res
                    break
            if lig_res:
                break

    if lig_res is None:
        return []

    # ── Ligand heavy-atom coordinates ─────────────────────────────────────────
    # Explicitly convert numpy.float32 → Python float so json.dump works.
    lig_coords: list[list[float]] = [
        [float(x) for x in atom.coord]
        for atom in lig_res.get_atoms()
        if atom.element and atom.element.upper() != "H"
    ]
    if not lig_coords:
        return []

    lig_centroid: list[float] = [
        float(sum(c[i] for c in lig_coords) / len(lig_coords))
        for i in range(3)
    ]

    # ── Process each residue token ─────────────────────────────────────────────
    records: list[PoseRecord] = []

    for token in entry.res_tokens:
        if not token or token[0].upper() not in _AA_CODES:
            continue

        res_code = token[0].upper()
        try:
            res_num = int(token[1:])
        except ValueError:
            continue

        # Locate residue — check receptor chain first, then any chain
        protein_res = None
        for chain_id in [entry.rec_chain] + [c.id for c in model if c.id != entry.rec_chain]:
            try:
                chain_obj = model[chain_id]
            except KeyError:
                continue
            for res in chain_obj:
                if res.id[1] == res_num and res.id[0] == " " and res.resname.strip():
                    # Verify AA 1-letter matches using local reverse map
                    aa_3 = res.resname.strip().upper()
                    aa_1 = _AA_3TO1.get(aa_3, "X")
                    if aa_1 == res_code:
                        protein_res = res
                        break
            if protein_res:
                break

        if protein_res is None:
            continue

        # Cα coordinate (float conversion for JSON serialization)
        ca_coord: list[float] | None = None
        if "CA" in protein_res:
            ca_coord = [float(x) for x in protein_res["CA"].coord]

        # Cα → ligand centroid distance
        ca_lig_dist: float | None = (
            float(_vec_dist(ca_coord, lig_centroid)) if ca_coord else None
        )

        # Minimum atom-to-atom distance (any protein residue atom ↔ any lig heavy atom)
        res_coords: list[list[float]] = [
            [float(x) for x in a.coord]
            for a in protein_res.get_atoms()
            if a.element and a.element.upper() != "H"
        ]
        if not res_coords:
            continue

        min_dist = min(
            _vec_dist(rc, lc)
            for rc in res_coords
            for lc in lig_coords
        )

        records.append(PoseRecord(
            pdb_id       = entry.pdb_id,
            lig_id       = entry.lig_id,
            lig_chain    = entry.lig_chain,
            lig_serial   = entry.lig_serial,
            rec_chain    = entry.rec_chain,
            res_code     = res_code,
            res_num      = res_num,
            fg_name      = fg_name,
            ca_coord     = ca_coord,
            lig_centroid = lig_centroid,
            ca_lig_dist  = round(ca_lig_dist, 3) if ca_lig_dist else None,
            min_atom_dist= round(min_dist, 3),
            resolution   = entry.resolution,
            ec_number    = entry.ec_number,
        ))

    return records


# ── Output 1: residue_3d_poses.json ──────────────────────────────────────────

class _SafeEncoder(json.JSONEncoder):
    """JSON encoder that converts numpy scalar types to Python natives."""
    def default(self, obj: object) -> object:
        try:
            import numpy as _np  # noqa: PLC0415
            if isinstance(obj, _np.floating):
                return float(obj)
            if isinstance(obj, _np.integer):
                return int(obj)
            if isinstance(obj, _np.ndarray):
                return obj.tolist()
        except ImportError:
            pass
        return super().default(obj)


def save_residue_3d_poses(all_records: dict[str, dict[str, list[PoseRecord]]]) -> None:
    """Save {fg: {residue: [PoseRecord, ...]}} to db/residue_3d_poses.json.

    Each PoseRecord is serialised as a plain dict (no dataclass-specific types).
    Uses _SafeEncoder to handle any residual numpy scalar types from BioPython.
    """
    out: dict = {}
    for fg, res_dict in all_records.items():
        out[fg] = {}
        for res, records in res_dict.items():
            out[fg][res] = [asdict(r) for r in records]

    with open(POSES_JSON, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, ensure_ascii=False, cls=_SafeEncoder)
    print(f"  Saved: {POSES_JSON}")


# ── Output 2: local_env SDF files ────────────────────────────────────────────

def _fg_safe(fg_name: str) -> str:
    """Convert FG name to a filesystem-safe string."""
    return re.sub(r"[^A-Za-z0-9]+", "_", fg_name).strip("_")


def save_local_env_sdf(
    all_records: dict[str, dict[str, list[PoseRecord]]],
    pdb_cache: dict[str, Path | None],
) -> None:
    """Write one SDF per FG × residue pair containing the ligand 3D structures.

    Each molecule in the SDF represents one binding event.
    Properties added per molecule:
      PDB_ID, LIG_ID, RES_CODE, RES_NUM, RESOLUTION,
      CA_LIG_DIST, MIN_ATOM_DIST, EC_NUMBER
    """
    try:
        from Bio.PDB import PDBParser
        from rdkit.Chem import AllChem
        from rdkit.Chem.rdmolfiles import SDWriter
    except ImportError as e:
        print(f"  ! SDF output requires BioPython + RDKit: {e}")
        return

    LOCAL_ENV_DIR.mkdir(parents=True, exist_ok=True)
    parser = PDBParser(QUIET=True)

    for fg, res_dict in all_records.items():
        fg_safe = _fg_safe(fg)
        for res, records in res_dict.items():
            sdf_path = LOCAL_ENV_DIR / f"{fg_safe}_{res}.sdf"
            writer = SDWriter(str(sdf_path))
            n_written = 0

            for rec in records:
                pdb_path = pdb_cache.get(rec.pdb_id)
                if pdb_path is None or not pdb_path.exists():
                    continue

                try:
                    structure = parser.get_structure("x", str(pdb_path))
                    model = structure[0]

                    # Find ligand residue (primary: lig_chain + lig_serial)
                    lig_res = None
                    for chain in model:
                        if chain.id != rec.lig_chain:
                            continue
                        for bres in chain:
                            if (bres.id[0].startswith("H_") and
                                    bres.id[1] == rec.lig_serial):
                                lig_res = bres
                                break
                        if lig_res:
                            break
                    # Fallback: first HETATM matching lig_id anywhere
                    if lig_res is None:
                        for chain in model:
                            for bres in chain:
                                if (bres.id[0].startswith("H_") and
                                        bres.resname.strip() == rec.lig_id.strip()):
                                    lig_res = bres
                                    break
                            if lig_res:
                                break

                    if lig_res is None:
                        continue

                    # Build RDKit mol from ligand heavy atoms + coords
                    rw_mol = Chem.RWMol()
                    conf_coords: list[tuple[float, float, float]] = []

                    for atom in lig_res.get_atoms():
                        if not atom.element or atom.element.upper() == "H":
                            continue
                        elem = atom.element.capitalize()
                        try:
                            rdatom = Chem.Atom(elem)
                            rw_mol.AddAtom(rdatom)
                            conf_coords.append(tuple(float(x) for x in atom.coord))
                        except Exception:
                            continue

                    if not conf_coords:
                        continue

                    from rdkit.Geometry import rdGeometry
                    conf = Chem.Conformer(len(conf_coords))
                    for i, (x, y, z) in enumerate(conf_coords):
                        conf.SetAtomPosition(i, rdGeometry.Point3D(x, y, z))
                    mol = rw_mol.GetMol()
                    mol.AddConformer(conf, assignId=True)

                    # Tag with metadata
                    mol.SetProp("PDB_ID",        rec.pdb_id)
                    mol.SetProp("LIG_ID",        rec.lig_id)
                    mol.SetProp("FG_NAME",       rec.fg_name)
                    mol.SetProp("RES_CODE",      rec.res_code)
                    mol.SetProp("RES_NUM",       str(rec.res_num))
                    mol.SetProp("RESOLUTION",    str(rec.resolution or "NMR"))
                    mol.SetProp("CA_LIG_DIST",   str(rec.ca_lig_dist or ""))
                    mol.SetProp("MIN_ATOM_DIST", str(rec.min_atom_dist))
                    mol.SetProp("EC_NUMBER",     rec.ec_number or "")

                    writer.write(mol)
                    n_written += 1

                except Exception:
                    continue

            writer.close()
            if n_written == 0 and sdf_path.exists():
                sdf_path.unlink()   # remove empty SDF
            elif n_written > 0:
                print(f"  SDF: {sdf_path.name}  ({n_written} mols)")


# ── Output 3: pharmacophore_stats.json ────────────────────────────────────────

def save_pharmacophore_stats(
    all_records: dict[str, dict[str, list[PoseRecord]]],
) -> None:
    """Compute distance statistics per FG × residue and save JSON.

    For each (fg, residue) pair, computes:
      - n_observations        total binding events with valid distance data
      - ca_lig_dist           Cα → ligand centroid (Å)  [mean, std, min, max, p25, p75]
      - min_atom_dist         min any-atom distance (Å) [same stats]
    """
    stats: dict = {}

    for fg, res_dict in all_records.items():
        stats[fg] = {}
        for res, records in res_dict.items():
            ca_dists   = [r.ca_lig_dist   for r in records if r.ca_lig_dist   is not None]
            min_dists  = [r.min_atom_dist for r in records if r.min_atom_dist is not None]

            def _desc(vals: list[float]) -> dict:
                if not vals:
                    return {"n": 0}
                arr = sorted(vals)
                n = len(arr)
                mean = sum(arr) / n
                std = math.sqrt(sum((v - mean) ** 2 for v in arr) / n) if n > 1 else 0.0
                p25 = arr[max(0, int(n * 0.25) - 1)]
                p75 = arr[min(n - 1, int(n * 0.75))]
                return {
                    "n":    n,
                    "mean": round(mean, 3),
                    "std":  round(std, 3),
                    "min":  round(arr[0], 3),
                    "max":  round(arr[-1], 3),
                    "p25":  round(p25, 3),
                    "p75":  round(p75, 3),
                }

            stats[fg][res] = {
                "n_observations": len(records),
                "ca_lig_dist":    _desc(ca_dists),
                "min_atom_dist":  _desc(min_dists),
            }

    with open(STATS_JSON, "w", encoding="utf-8") as f:
        json.dump(stats, f, indent=2, ensure_ascii=False, cls=_SafeEncoder)
    print(f"  Saved: {STATS_JSON}")


# ── Main pipeline ──────────────────────────────────────────────────────────────

def run(
    top_pairs:      int        = DEFAULT_TOP_PAIRS,
    top_per_fg:     int | None = None,
    max_per_pair:   int        = DEFAULT_MAX_PER_PAIR,
    max_resolution: float      = DEFAULT_MAX_RESOL,
    fg_filter:      str | None = None,
    res_filter:     str | None = None,
    dry_run:        bool       = False,
    keep_pdb:       bool       = False,
) -> None:
    """Full pose extraction pipeline."""

    # ── Step 1: Parse BioLiP ──────────────────────────────────────────────────
    if not BIOLIP_NR_CACHE.exists():
        print(f"ERROR: BioLiP cache not found at {BIOLIP_NR_CACHE}")
        print("  Run: python utils/interaction_analyzer.py --local db/BioLiP_nr.txt.gz")
        sys.exit(1)

    print("Parsing BioLiP (full metadata) ...")
    entries = parse_biolip_full(BIOLIP_NR_CACHE)
    print(f"  {len(entries):,} entries parsed")

    # ── Step 2: Load SMILES cache & build FG × residue index ─────────────────
    print("Loading SMILES cache ...")
    smiles_cache = _load_smiles_cache()
    print(f"  {len(smiles_cache):,} ligand SMILES loaded")

    print("Building FG × residue index ...")
    fg_res_index = build_fg_residue_index(entries, smiles_cache)

    # ── Step 3: Select representatives ────────────────────────────────────────
    print(f"\nSelecting representatives "
          f"(top {top_pairs} pairs, max {max_per_pair}/pair, "
          f"res <= {max_resolution} A) ...")
    selected = select_representatives(
        fg_res_index, top_pairs, max_per_pair, max_resolution,
        fg_filter, res_filter,
        top_per_fg=top_per_fg,
    )

    total_structs = sum(
        len(ents)
        for res_dict in selected.values()
        for ents in res_dict.values()
    )
    unique_pdbs = {
        e.pdb_id
        for res_dict in selected.values()
        for ents in res_dict.values()
        for e in ents
    }
    print(f"\n  Total records selected : {total_structs:,}")
    print(f"  Unique PDB structures  : {len(unique_pdbs):,}")

    if dry_run:
        print("\n[dry-run] Stopping before downloads.")
        return

    # ── Step 4: Download PDB structures ───────────────────────────────────────
    print(f"\nDownloading {len(unique_pdbs)} PDB structures to {TMP_DIR}/ ...")
    pdb_cache: dict[str, Path | None] = {}
    n_ok = n_fail = 0

    for i, pdb_id in enumerate(sorted(unique_pdbs), 1):
        path = download_pdb(pdb_id)
        pdb_cache[pdb_id] = path
        if path:
            n_ok += 1
        else:
            n_fail += 1
        if i % 50 == 0:
            print(f"  {i}/{len(unique_pdbs)} downloaded ({n_fail} failed) ...")

    print(f"  Downloads: {n_ok} OK, {n_fail} failed")

    # ── Step 5: Extract 3D pose data ──────────────────────────────────────────
    print("\nExtracting 3D pose data ...")
    all_records: dict[str, dict[str, list[PoseRecord]]] = defaultdict(
        lambda: defaultdict(list)
    )
    n_extracted = n_empty = 0

    for fg, res_dict in selected.items():
        for res, ents in res_dict.items():
            for entry in ents:
                pdb_path = pdb_cache.get(entry.pdb_id)
                if pdb_path is None or not pdb_path.exists():
                    continue
                records = extract_pose_records(pdb_path, entry, fg)
                res_records = [r for r in records if r.res_code == res]
                if res_records:
                    all_records[fg][res].extend(res_records)
                    n_extracted += len(res_records)
                else:
                    n_empty += 1

    print(f"  Extracted: {n_extracted:,} pose records "
          f"({n_empty} entries yielded no matching residue)")

    # ── Step 6: Save outputs ──────────────────────────────────────────────────
    print("\nSaving outputs ...")
    save_residue_3d_poses(all_records)
    save_local_env_sdf(all_records, pdb_cache)
    save_pharmacophore_stats(all_records)

    # ── Step 7: Cleanup ───────────────────────────────────────────────────────
    if not keep_pdb and TMP_DIR.exists():
        n_deleted = 0
        for pdb_file in TMP_DIR.glob("*.pdb"):
            pdb_file.unlink()
            n_deleted += 1
        print(f"\n  Deleted {n_deleted} PDB files from {TMP_DIR}/")

    print("\nDone.")


# ── CLI ────────────────────────────────────────────────────────────────────────

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Extract 3D residue poses from BioLiP representative structures."
    )
    p.add_argument(
        "--top-pairs", type=int, default=DEFAULT_TOP_PAIRS, metavar="N",
        help=f"Process top N FG x residue pairs globally (default: {DEFAULT_TOP_PAIRS}). "
             "Ignored if --top-per-fg is set.",
    )
    p.add_argument(
        "--top-per-fg", type=int, default=None, metavar="N",
        help="Top N residue partners PER FG (recommended for full coverage). "
             "Example: --top-per-fg 3 -> 21 FGs x 3 residues = up to 63 pairs.",
    )
    p.add_argument(
        "--max-per-pair", type=int, default=DEFAULT_MAX_PER_PAIR, metavar="N",
        help=f"Max representative structures per pair (default: {DEFAULT_MAX_PER_PAIR})",
    )
    p.add_argument(
        "--resolution", type=float, default=DEFAULT_MAX_RESOL, metavar="A",
        help=f"Max resolution to accept in Angstrom (default: {DEFAULT_MAX_RESOL}). "
             "NMR structures always accepted.",
    )
    p.add_argument(
        "--fg", default=None, metavar="FG_NAME",
        help='Filter to a single FG, e.g. "Phenol" (default: all FGs)',
    )
    p.add_argument(
        "--residue", default=None, metavar="AA",
        help="Filter to a single residue 1-letter code, e.g. Y (default: all)",
    )
    p.add_argument(
        "--dry-run", action="store_true",
        help="Show selection plan without downloading structures",
    )
    p.add_argument(
        "--keep-pdb", action="store_true",
        help="Keep downloaded PDB files in tmp/ after extraction",
    )
    return p


def main() -> None:
    """CLI entry point."""
    args = _build_parser().parse_args()
    run(
        top_pairs      = args.top_pairs,
        top_per_fg     = args.top_per_fg,
        max_per_pair   = args.max_per_pair,
        max_resolution = args.resolution,
        fg_filter      = args.fg,
        res_filter     = args.residue,
        dry_run        = args.dry_run,
        keep_pdb       = args.keep_pdb,
    )


if __name__ == "__main__":
    main()
