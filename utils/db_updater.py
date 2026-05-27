"""Build and refresh db/fg_database.json.

Source priority per field:
  description  → ChEBI (ontology-level definition) → PubChem fallback
  pubchem_cid  → PubChem name search
  chebi_id     → ChEBI known-ID table → ChEBI name search

Usage:
    python utils/db_updater.py
    python utils/db_updater.py --fg-only "Carboxylic acid (COOH),Ester"
    python utils/db_updater.py --dry-run

APIs used:
  PubChem: https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest
  ChEBI:   https://www.ebi.ac.uk/chebi/searchId.do (HTML, parsed)
"""

import argparse
import json
import re
import sys
import time
from datetime import date
from pathlib import Path

import requests

# ── Paths ─────────────────────────────────────────────────────────────────────
_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT))

from constants.fg_smarts import FG_SMARTS  # noqa: E402

DB_PATH = _ROOT / "db" / "fg_database.json"

# ── Rate limits ────────────────────────────────────────────────────────────────
PUBCHEM_DELAY = 0.25  # PubChem: max 5 req/sec for anonymous
CHEBI_DELAY   = 0.5   # ChEBI: be polite

# ── ChEBI: known IDs for our functional group set ─────────────────────────────
# These are stable ChEBI ontology IDs; unlikely to change.
_CHEBI_IDS: dict[str, str] = {
    "Carboxylic acid (COOH)": "CHEBI:33575",   # carboxylic acid
    "Ester":                  "CHEBI:35701",   # ester
    "Ether (C-O-C)":          "CHEBI:25698",   # ether (ROR)
    "Aromatic ring":          "CHEBI:22712",   # benzenes (aromatic ring class)
    "Amide (CONH)":           "CHEBI:37622",   # carboxamide
    "Amine (NH2)":            "CHEBI:32952",   # amine
    "Carbonyl (C=O)":         "CHEBI:23019",   # carbonyl group (no def, but correct)
    "Hydroxyl (OH)":          "CHEBI:43176",   # hydroxy group (no def, but correct)
    "Halogen":                "CHEBI:17792",   # organohalogen compound
    "Nitrile (CN)":           "CHEBI:18379",   # nitrile RC≡N
    "Nitro (NO2)":            "CHEBI:35716",   # C-nitro compound
    "Sulfonamide":            "CHEBI:35358",   # sulfonamide
    "Imidazole":              "CHEBI:24780",   # imidazoles class
    "Benzene":                "CHEBI:16716",   # benzene
    "Tertiary amine":         "CHEBI:32876",   # tertiary amine
}

# ── Manual descriptions for FGs with no usable API definition ─────────────────
# Used when ChEBI has no definition and PubChem returns an irrelevant compound.
_MANUAL_DESCRIPTIONS: dict[str, str] = {
    "Carbonyl (C=O)": (
        "A functional group consisting of a carbon atom doubly bonded to an oxygen atom "
        "(C=O). Found in aldehydes, ketones, carboxylic acids, esters, and amides."
    ),
    "Hydroxyl (OH)": (
        "A functional group consisting of an oxygen atom bonded to a hydrogen atom (‒OH). "
        "Present in alcohols and phenols; confers hydrogen-bond donor capability."
    ),
}

# ── Seed: known target classes ─────────────────────────────────────────────────
# Cannot be derived automatically from any API.
# Sourced from pharmacology / medicinal chemistry literature.
_TARGET_CLASS_SEEDS: dict[str, list[str]] = {
    "Carboxylic acid (COOH)": ["protease", "transporter", "nuclear receptor"],
    "Ester":                  ["esterase", "lipase"],
    "Ether (C-O-C)":          [],
    "Aromatic ring":          ["kinase", "GPCR", "nuclear receptor"],
    "Amide (CONH)":           ["protease", "kinase"],
    "Amine (NH2)":            ["GPCR", "ion channel", "transporter"],
    "Carbonyl (C=O)":         ["oxidoreductase", "aldehyde dehydrogenase"],
    "Hydroxyl (OH)":          ["kinase", "phosphatase", "transporter"],
    "Halogen":                ["kinase", "ion channel"],
    "Nitrile (CN)":           ["protease", "nitrile hydratase"],
    "Nitro (NO2)":            ["nitroreductase"],
    "Sulfonamide":            ["carbonic anhydrase", "COX"],
    "Imidazole":              ["histamine receptor", "cytochrome P450"],
    "Benzene":                ["kinase", "GPCR", "nuclear receptor"],
    "Tertiary amine":         ["GPCR", "ion channel", "cholinesterase"],
}

# ── ChEBI helpers ──────────────────────────────────────────────────────────────

_CHEBI_DEF_RE = re.compile(
    r'Definition</td><td[^>]*>(.*?)</td>', re.DOTALL
)
_CHEBI_ID_SEARCH_RE = re.compile(r'CHEBI:\d+')


def _chebi_fetch_by_id(chebi_id: str) -> dict:
    """Fetch definition and confirm ID from a known ChEBI ID.

    Parses the ChEBI entity HTML page.
    Returns dict with 'chebi_id' and 'description' (either may be None).
    """
    url = f"https://www.ebi.ac.uk/chebi/searchId.do?chebiId={chebi_id}"
    try:
        r = requests.get(url, timeout=15)
        time.sleep(CHEBI_DELAY)
        if r.status_code != 200:
            return {}
        m = _CHEBI_DEF_RE.search(r.text)
        if m:
            raw = m.group(1).strip()
            # Strip any HTML tags inside the definition (rare but possible)
            definition = re.sub(r'<[^>]+>', '', raw).strip()
            return {"chebi_id": chebi_id, "description": definition}
        # Page loaded but no definition block found
        return {"chebi_id": chebi_id, "description": None}
    except Exception:
        return {}


def _chebi_search_by_name(fg_name: str) -> dict:
    """Search ChEBI by name and extract the first ChEBI ID from results HTML.

    Used as fallback when a known ID is not in _CHEBI_IDS.
    Returns dict with 'chebi_id' and optionally 'description'.
    """
    clean = re.sub(r'\s*\(.*?\)', '', fg_name).strip()
    url = f"https://www.ebi.ac.uk/chebi/search?q={requests.utils.quote(clean)}"
    try:
        r = requests.get(url, timeout=15)
        time.sleep(CHEBI_DELAY)
        if r.status_code != 200:
            return {}
        ids = _CHEBI_ID_SEARCH_RE.findall(r.text)
        if not ids:
            return {}
        first_id = ids[0]
        return _chebi_fetch_by_id(first_id)
    except Exception:
        return {}


def _fetch_chebi(fg_name: str) -> dict:
    """Fetch ChEBI data for a functional group.

    Uses pre-seeded ID if available; falls back to name search.
    """
    chebi_id = _CHEBI_IDS.get(fg_name)
    if chebi_id:
        return _chebi_fetch_by_id(chebi_id)
    return _chebi_search_by_name(fg_name)


# ── PubChem helpers ────────────────────────────────────────────────────────────

def _pubchem_name_variants(fg_name: str) -> list[str]:
    """Return name variants to try in PubChem, best-first.

    Strips parentheticals (e.g. 'Ether (C-O-C)' → 'Ether') and tries
    first-word only as a last resort.
    """
    stripped = re.sub(r'\s*\(.*?\)', '', fg_name).strip()
    candidates = [fg_name]
    if stripped and stripped != fg_name:
        candidates.append(stripped)
    first = stripped.split()[0] if stripped else ""
    if first and first not in candidates and len(first) > 3:
        candidates.append(first)
    return candidates


def _pubchem_cid_by_name(name: str) -> int | None:
    """Return the first PubChem CID matching this name, or None."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{requests.utils.quote(name)}/cids/JSON"
    try:
        r = requests.get(url, timeout=10)
        time.sleep(PUBCHEM_DELAY)
        if r.status_code == 200:
            return r.json()["IdentifierList"]["CID"][0]
    except Exception:
        pass
    return None


def _pubchem_description(cid: int) -> str | None:
    """Return the first non-empty description for a PubChem CID, or None."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/description/JSON"
    try:
        r = requests.get(url, timeout=10)
        time.sleep(PUBCHEM_DELAY)
        if r.status_code == 200:
            for section in r.json().get("InformationList", {}).get("Information", []):
                desc = section.get("Description", "").strip()
                if desc:
                    return desc
    except Exception:
        pass
    return None


def _fetch_pubchem(fg_name: str) -> dict:
    """Fetch PubChem CID (and description if available) for a functional group."""
    result: dict = {"pubchem_cid": None}

    cid: int | None = None
    for candidate in _pubchem_name_variants(fg_name):
        cid = _pubchem_cid_by_name(candidate)
        if cid is not None:
            break

    if cid is None:
        return result

    result["pubchem_cid"] = cid

    desc = _pubchem_description(cid)
    if desc:
        result["pubchem_description"] = desc  # kept separate; ChEBI takes priority

    return result


# ── Core build logic ───────────────────────────────────────────────────────────

def build_entry(fg_name: str, smarts: str) -> dict:
    """Build a complete fg_database.json entry for one functional group.

    Field priority:
      description → ChEBI first, PubChem fallback
      pubchem_cid → PubChem
      chebi_id    → ChEBI
    """
    entry: dict = {
        "smarts":               smarts,
        "pubchem_cid":          None,
        "chebi_id":             None,
        "description":          None,
        "known_target_classes": _TARGET_CLASS_SEEDS.get(fg_name, []),
        "source":               None,
    }

    # 1. ChEBI — ontology-level description (best quality)
    chebi_data = _fetch_chebi(fg_name)
    if chebi_data.get("chebi_id"):
        entry["chebi_id"] = chebi_data["chebi_id"]
    if chebi_data.get("description"):
        entry["description"] = chebi_data["description"]
        entry["source"] = "chebi"

    # 2. Manual curated — overrides PubChem for FGs with known-bad API results
    if fg_name in _MANUAL_DESCRIPTIONS and entry["description"] is None:
        entry["description"] = _MANUAL_DESCRIPTIONS[fg_name]
        entry["source"] = "manual"

    # 3. PubChem — compound CID (always); description only if nothing better found
    pubchem_data = _fetch_pubchem(fg_name)
    if pubchem_data.get("pubchem_cid"):
        entry["pubchem_cid"] = pubchem_data["pubchem_cid"]
    if entry["description"] is None and pubchem_data.get("pubchem_description"):
        entry["description"] = pubchem_data["pubchem_description"]
        entry["source"] = "pubchem"

    return entry


def build_database(target_fgs: list[str] | None = None) -> dict:
    """Build the complete fg_database.json structure.

    Args:
        target_fgs: Limit update to these FG names. None = update all.

    Returns:
        Complete database dict ready for JSON serialisation.
    """
    fg_subset = {
        name: smarts
        for name, smarts in FG_SMARTS.items()
        if target_fgs is None or name in target_fgs
    }

    # Preserve existing entries that are not being updated
    existing: dict = {}
    if DB_PATH.exists():
        with open(DB_PATH, encoding="utf-8") as f:
            existing = json.load(f)

    functional_groups: dict = existing.get("functional_groups", {})

    total = len(fg_subset)
    for i, (fg_name, smarts) in enumerate(fg_subset.items(), start=1):
        print(f"  [{i}/{total}] {fg_name} ...", end=" ", flush=True)
        entry = build_entry(fg_name, smarts)
        functional_groups[fg_name] = entry
        print(entry["source"] or "not found")

    return {
        "last_updated":      str(date.today()),
        "sources":           ["chebi", "pubchem"],
        "functional_groups": functional_groups,
    }


def save_database(db: dict) -> None:
    """Write the database to db/fg_database.json."""
    DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(DB_PATH, "w", encoding="utf-8") as f:
        json.dump(db, f, indent=2, ensure_ascii=False)
    print(f"\nSaved: {DB_PATH}")


# ── Entry point ────────────────────────────────────────────────────────────────

def main() -> None:
    """Parse args and run the database update."""
    # Windows consoles default to cp950/cp932; force UTF-8 for ChEBI descriptions.
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

    parser = argparse.ArgumentParser(description="Update fg_database.json")
    parser.add_argument(
        "--fg-only",
        default=None,
        help='Comma-separated FG names to update, e.g. "Ester,Halogen"',
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print result without saving",
    )
    args = parser.parse_args()

    target_fgs = (
        [s.strip() for s in args.fg_only.split(",")]
        if args.fg_only
        else None
    )

    print("Building fg_database.json ...")
    db = build_database(target_fgs=target_fgs)

    if args.dry_run:
        print("\n-- DRY RUN (not saved) --")
        print(json.dumps(db, indent=2, ensure_ascii=False))
    else:
        save_database(db)
        print("Done.")


if __name__ == "__main__":
    main()
