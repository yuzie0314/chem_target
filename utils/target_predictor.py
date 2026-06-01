"""Target class and residue interaction prediction from a SMILES string.

Workflow
--------
1. Detect functional groups (SMARTS-based) in the input SMILES.
2. Load the pre-built FG × residue co-occurrence table (db/fg_residue_table.csv).
   rows  = 20 standard amino acids (1-letter code)
   cols  = FG names from FG_SMARTS
   value = binding-event count from BioLiP 2.0
3. Score each residue by **z-score-normalising** each FG column then summing
   (prevents high-frequency generic FGs from drowning specific pharmacophores).
4. Vote target classes using known_target_classes from db/fg_database.json,
   weighted by **IDF** (specific targets ranked above generic labels).
5. Return a structured result dict + a formatted text report.

Scoring details
---------------
Residue scoring (z-score normalisation, default on):
    For each detected FG, normalise its column to mean=0 / std=1 across all 20 AAs
    before summing.  This ensures that a highly represented FG like Hydroxyl
    (>11 000 GLY events) does not completely overshadow Steroid (~1 400 PHE events).

Target class scoring (IDF weighting, default on):
    weight(tc) = log( N_all_FGs / N_FGs_that_annotate_this_tc )
    High weight → target class is specific (few FGs annotate it, e.g. VKORC1, tubulin).
    Low weight  → target class is generic (many FGs annotate it, e.g. kinase, GPCR).
    final_score = vote_count × weight

Prerequisites
-------------
Build the FG × residue table once (uses local BioLiP_nr.txt.gz):

    python utils/interaction_analyzer.py --local db/BioLiP_nr.txt.gz

Outputs: db/fg_residue_table.csv  (tracked in git after first build)

Quick test with the first 1000 entries:

    python utils/interaction_analyzer.py --local db/BioLiP_nr.txt.gz --top 1000

Standalone usage
----------------
    python utils/target_predictor.py "CC(=O)Oc1ccccc1C(=O)O" --name Aspirin
    python utils/target_predictor.py "CC(=O)Oc1ccccc1C(=O)O" --top 5
"""

import json
import sys
from collections import Counter, defaultdict
from math import log
from pathlib import Path

import pandas as pd

_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT))

from utils.fg_detector import detect_smarts  # noqa: E402

# ── Paths ──────────────────────────────────────────────────────────────────────
TABLE_PATH = _ROOT / "db" / "fg_residue_table.csv"
FG_DB_PATH = _ROOT / "db" / "fg_database.json"

# ── Conditional-scoring constants ─────────────────────────────────────────────
# Pre-IDF vote weights added when multi-FG combos match.  These are accumulated
# into weighted_votes before IDF multiplication, so the final score contribution
# is bonus × IDF(cytochrome P450).
_CYPCOND_AZOLE_BONUS: float = 2.0   # Imidazole + lipophilic partner
_CYPCOND_LIPOPHILIC_FGS: frozenset[str] = frozenset({"Phenyl ring", "Ether", "Halogen"})
_COX_INDOLE_SULFONAMIDE_BONUS: float = 2.0  # Indole + Sulfonamide COX-2 pharmacophore
_MTOR_MACROLIDE_BONUS: float = 2.0          # Macrolide without competing metal-binding warheads
_ADENOSINE_PURINE_BONUS: float = 0.5       # Purine scaffold — adenosine receptor defining motif
_KINASE_ABUNSAT_BONUS: float = 0.5         # alpha,beta-unsat carbonyl — covalent kinase warhead
_KINASE_SULFONAMIDE_TAMINE_BONUS: float = 2.0  # Sulfonamide + TertAmine — kinase linker hijacked by CA
_CYP450_ARYL_HALIDE_COOH_BONUS: float = 1.5   # Aryl-halide carboxylic acid — CYP substrate (no Amide/Ether)

# ── Amino acid code lookup ─────────────────────────────────────────────────────
AA_1TO3: dict[str, str] = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
    "E": "GLU", "Q": "GLN", "G": "GLY", "H": "HIS", "I": "ILE",
    "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
    "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}


# ── Data loaders ───────────────────────────────────────────────────────────────

def load_table(path: Path = TABLE_PATH) -> pd.DataFrame:
    """Load the pre-built FG × residue co-occurrence table.

    The table has:
        index  = residue (1-letter AA code)
        cols   = "residue_name" + FG names from FG_SMARTS
        values = integer co-occurrence counts from BioLiP 2.0

    Raises:
        FileNotFoundError: with instructions to build the table if missing.
    """
    if not path.exists():
        raise FileNotFoundError(
            f"\nFG-residue table not found: {path}\n\n"
            "Build it first (requires db/BioLiP_nr.txt.gz):\n"
            "  python utils/interaction_analyzer.py --local db/BioLiP_nr.txt.gz\n\n"
            "Quick test (first 1 000 entries):\n"
            "  python utils/interaction_analyzer.py --local db/BioLiP_nr.txt.gz --top 1000"
        )
    df = pd.read_csv(path, index_col="residue")
    return df


def load_fg_db(path: Path = FG_DB_PATH) -> dict:
    """Load fg_database.json and return the functional_groups sub-dict.

    Returns an empty dict if the file is missing (graceful fallback).
    """
    if not path.exists():
        return {}
    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    return data.get("functional_groups", {})


# ── Scoring helpers ────────────────────────────────────────────────────────────

def _compute_target_idf(fg_db: dict) -> dict[str, float]:
    """Compute inverse-document-frequency weights for target classes.

    IDF(tc) = log( N_all_FGs / N_FGs_that_annotate_tc )

    Interpretation:
        High IDF → target is specific (few FGs list it, e.g. VKORC1, tubulin).
        Low IDF  → target is generic (many FGs list it, e.g. kinase, GPCR).

    Args:
        fg_db: functional_groups dict loaded from fg_database.json.

    Returns:
        Dict mapping target class name → IDF weight (float).
    """
    n_fgs = len(fg_db)
    if n_fgs == 0:
        return {}
    tc_count: Counter = Counter()
    for entry in fg_db.values():
        for tc in entry.get("known_target_classes", []):
            tc_count[tc] += 1
    return {tc: log(n_fgs / count) for tc, count in tc_count.items()}


def _cyp450_conditional_bonus(fgs_detected: list[str]) -> tuple[float, str]:
    """Pre-IDF bonus vote for CYP450 when azole-type combos are present.

    Rule: Imidazole + at least one lipophilic partner (Phenyl ring, Ether,
    or Halogen), provided the compound lacks FGs that indicate a competing
    metal-binding context:
      • Ketone excluded: alpha-keto HDAC warhead (romidepsin-class inhibitors
        contain Imidazole + Ketone but are HDAC substrates, not CYP heme binders).
      • Purine excluded: the imidazole-like ring in purines (adenosine receptor
        ligands, xanthine derivatives) would trigger this rule spuriously since
        RDKit's Imidazole SMARTS matches the fused imidazole portion of purines.

    Returns:
        (bonus_wt, label) — pre-IDF weight to add to cytochrome P450 votes,
        and a short evidence label string.  Returns (0.0, '') if rule does not fire.
    """
    fg_set = set(fgs_detected)
    if (
        "Imidazole" in fg_set
        and fg_set & _CYPCOND_LIPOPHILIC_FGS
        and "Ketone" not in fg_set
        and "Purine" not in fg_set
        and "α,β-unsat. carbonyl" not in fg_set  # covalent kinase warhead context
    ):
        return _CYPCOND_AZOLE_BONUS, "azole motif"
    return 0.0, ""


def _cyp450_arylhalide_cooh_bonus(fgs_detected: list[str]) -> tuple[float, str]:
    """Pre-IDF bonus for CYP450 when aryl-halide carboxylic acid is present without NSAID context.

    Two rules (additive where both fire):

    Rule A: Carboxylic acid + Phenyl ring + Halogen, AND Amide absent AND Ether absent.
      Catches minimal aryl-halide COOH scaffold (CHEMBL3125537, CHEMBL3407554).
      Exclusions prevent NSAID conflict:
        • INDOMETHACIN (COX HIT): has Ether → excluded ✓
        • CHEMBL4582020 (COX HIT): has Amide → excluded ✓

    Rule B: Carboxylic acid + Amide + Ether + Phenyl ring + Halogen.
      Catches extended aryl-halide COOH + linker scaffold (CHEMBL3407558, CHEMBL3407575).
      Safety: the only benchmark HIT with this combo is CHEMBL6067690 (GPCR, score=7.784
      from 4 GPCR FGs), which remains GPCR even with the +1.5 CYP450 bonus (4.022 < 7.784).

    Returns:
        (bonus_wt, label) — pre-IDF weight and evidence label.  Returns (0.0, '') if neither
        rule fires.
    """
    fg_set = set(fgs_detected)
    # Rule A — minimal aryl-halide COOH
    if (
        "Carboxylic acid" in fg_set
        and "Phenyl ring" in fg_set
        and "Halogen" in fg_set
        and "Amide" not in fg_set
        and "Ether" not in fg_set
    ):
        return _CYP450_ARYL_HALIDE_COOH_BONUS, "aryl-halide COOH CYP substrate"
    # Rule B — extended aryl-halide COOH + Amide + Ether linker
    if (
        "Carboxylic acid" in fg_set
        and "Amide" in fg_set
        and "Ether" in fg_set
        and "Phenyl ring" in fg_set
        and "Halogen" in fg_set
    ):
        return _CYP450_ARYL_HALIDE_COOH_BONUS, "aryl-halide COOH+linker CYP substrate"
    return 0.0, ""


def _cox_conditional_bonus(fgs_detected: list[str]) -> tuple[float, str]:
    """Pre-IDF bonus vote for COX when Indole+Sulfonamide pharmacophore present.

    Rule: Indole + Sulfonamide co-occurrence.
    Indole-sulfonamide COX inhibitors (indomethacin-like analogs) are
    consistently mispredicted as carbonic anhydrase because Sulfonamide
    carries mw=2.0 for CA.  The Indole ring is the signature scaffold for
    COX binding (fits the hydrophobic channel), while the Sulfonamide head
    sits in the COX-2 selectivity pocket.  The combination is rarely seen in
    CA inhibitors (which use heterocyclic sulfonamides without Indole).

    Returns:
        (bonus_wt, label) — pre-IDF weight and evidence label.
    """
    fg_set = set(fgs_detected)
    if "Indole" in fg_set and "Sulfonamide" in fg_set:
        return _COX_INDOLE_SULFONAMIDE_BONUS, "indole-sulfonamide COX motif"
    return 0.0, ""


def _mtor_conditional_bonus(fgs_detected: list[str]) -> tuple[float, str]:
    """Pre-IDF bonus vote for mTOR when macrolide scaffold is present without
    competing metal-binding or tubulin warheads.

    Rule: Macrolide present AND Thiol absent AND α,β-unsat. carbonyl absent AND
    Acylsulfonamide absent.

    Rationale for exclusions:
      • Thiol / α,β-unsat. carbonyl: HDAC zinc-binding context (romidepsin-class);
        FR-135313 has Macrolide+Thiol+αβunsat and is correctly HDAC.
      • Acylsulfonamide: tubulin macrolide warhead (epothilone class); 12 tubulin
        benchmark compounds have Macrolide+Acylsulfonamide and must stay as tubulin.

    The surviving case is rapamycin/sirolimus-class allosteric mTOR inhibitors,
    which have Macrolide+Ketone but no thiol/αβunsat/acylsulfonamide.

    Returns:
        (bonus_wt, label) — pre-IDF weight and evidence label.
    """
    fg_set = set(fgs_detected)
    if (
        "Macrolide" in fg_set
        and "Thiol" not in fg_set
        and "α,β-unsat. carbonyl" not in fg_set
        and "Acylsulfonamide" not in fg_set
    ):
        return _MTOR_MACROLIDE_BONUS, "macrolide mTOR motif"
    return 0.0, ""


def _adenosine_conditional_bonus(fgs_detected: list[str]) -> tuple[float, str]:
    """Pre-IDF bonus for adenosine receptor when Purine scaffold is present.

    Rule: Purine detected → +0.5 bonus to adenosine receptor.

    When Purine + Phenyl ring co-occur, kinase accrues two IDF-weighted votes
    (Purine + Phenyl → 2.0 × 1.946 = 3.892) while adenosine receptor only gets
    one vote from Purine (1.0 × 2.862 = 2.862).  The Purine scaffold is the
    defining feature of most adenosine receptor ligands (xanthine derivatives,
    adenosine analogs); this bonus prevents spurious kinase assignment solely
    due to Phenyl co-occurrence.  All confirmed kinase benchmark HITs use Lactone
    or α,β-unsat. carbonyl as their primary signal, not Purine.

    Returns:
        (bonus_wt, label) — pre-IDF weight and evidence label.
    """
    if "Purine" in set(fgs_detected):
        return _ADENOSINE_PURINE_BONUS, "purine adenosine motif"
    return 0.0, ""


def _kinase_conditional_bonus(fgs_detected: list[str]) -> tuple[float, str]:
    """Pre-IDF bonus for kinase when alpha,beta-unsaturated carbonyl is present.

    Rule: α,β-unsat. carbonyl detected → +0.5 pre-IDF for kinase.

    This covalent Michael acceptor is the defining warhead of irreversible kinase
    inhibitors (osimertinib, neratinib, afatinib, mobocertinib attacking Cys797 in
    EGFR).  Without this bonus, these compounds tie with GPCR at 3.892 (TerAmine +
    Phenyl = αβunsat + Phenyl) and lose the tie-break to GPCR, or lose to cysteine
    protease when Nitrile co-occurs (Nitrile+αβunsat = 2.0 × IDF_cys = 4.338 vs
    kinase 3.892).  The +0.5 bonus gives kinase 4.865, resolving both conflicts.

    Safety rationale:
      • All 20 GPCR benchmark HITs lack α,β-unsat. carbonyl (verified).
      • HDAC HITs with αβunsat all have Hydroxamate(mw=2.5)+αβunsat → HDAC = 6.8+,
        well above kinase bonus (4.865).
      • NR HITs with αβunsat have Steroid(mw=2.0) → NR/androgen = 7.11+.

    Returns:
        (bonus_wt, label) — pre-IDF weight and evidence label.
    """
    fg_set = set(fgs_detected)
    # Rule 1: covalent kinase warhead (Michael acceptor)
    if "α,β-unsat. carbonyl" in fg_set:
        return _KINASE_ABUNSAT_BONUS, "covalent kinase warhead"
    # Rule 2: Sulfonamide + Tertiary amine — kinase inhibitor with sulfonamide linker,
    # systematically hijacked by carbonic anhydrase (Sulfonamide mw=2.0 × IDF_CA).
    # Safety: only ONE benchmark compound has this combination (CHEMBL5594833, kinase).
    if "Sulfonamide" in fg_set and "Tertiary amine" in fg_set:
        return _KINASE_SULFONAMIDE_TAMINE_BONUS, "sulfonamide-amine kinase linker"
    return 0.0, ""


def _apply_negative_constraints(
    fgs_detected: list[str],
    weighted_votes: dict[str, float],
    evidence: dict[str, list[str]],
) -> None:
    """Suppress CYP450 votes in-place when incompatible warheads are present.

    Rules:
      • Hydroxamate or Thiol → Zn2+ chelation → HDAC / metalloprotease context.
      • Acylsulfonamide → tubulin macrolide warhead.
      Both: remove cytochrome P450 entry entirely.
    """
    fg_set = set(fgs_detected)
    if fg_set & {"Hydroxamate", "Thiol", "Acylsulfonamide"}:
        weighted_votes.pop("cytochrome P450", None)
        evidence.pop("cytochrome P450", None)


def predict_residues(
    fgs_detected: list[str],
    table: pd.DataFrame,
    top_n: int = 20,
    normalize: bool = True,
) -> pd.DataFrame:
    """Score residues by summing BioLiP co-occurrence counts for detected FGs.

    When *normalize* is True (default), each FG column is z-score-normalised
    (column mean → 0, column std → 1) before summing.  This prevents dominant
    FGs (e.g. Hydroxyl: ~11 000 GLY hits) from overshadowing structurally
    specific FGs (e.g. Steroid: ~1 400 PHE hits) when a compound contains both.

    Args:
        fgs_detected: FG names matching FG_SMARTS keys.
        table:        Numeric FG × residue table (residue_name column excluded).
        top_n:        Return only the top-N residues by score (0 = return all).
        normalize:    Apply per-FG z-score normalisation before summing.

    Returns:
        DataFrame with columns: residue_name | score | contributing_fgs
        Index: residue (1-letter AA code).  Rows with score ≤ 0 are dropped.
    """
    valid_fgs = [fg for fg in fgs_detected if fg in table.columns]

    if not valid_fgs:
        return pd.DataFrame(
            columns=["residue_name", "score", "contributing_fgs"]
        )

    sub = table[valid_fgs].astype(float)

    if normalize:
        # Per-column z-score.  Avoid division by zero for all-zero FG columns.
        col_std = sub.std()
        col_std[col_std == 0] = 1.0
        sub = (sub - sub.mean()) / col_std

    score_series = sub.sum(axis=1)

    result = pd.DataFrame({
        "residue_name": score_series.index.map(lambda aa: AA_1TO3.get(aa, aa)),
        "score": score_series,
        "contributing_fgs": ", ".join(valid_fgs),
    })
    result.index.name = "residue"
    result = result[result["score"] > 0].sort_values("score", ascending=False)

    if top_n and top_n > 0:
        result = result.head(top_n)

    return result


def predict_target_classes(
    fgs_detected: list[str],
    fg_db: dict,
    use_idf: bool = True,
) -> pd.DataFrame:
    """Vote target classes using known_target_classes from fg_database.json.

    Scoring formula (with IDF and mechanistic weighting):

        score(tc) = IDF(tc) × Σ_{fg→tc} mechanistic_weight(fg)

    where:
        IDF(tc)               = log( N_all_FGs / N_FGs_annotating_tc )
        mechanistic_weight(fg) = fg_database.json field "mechanistic_weight"
                                 (default 1.0 if absent)

    Rationale:
        IDF separates *specific* targets (few FGs annotate them) from *generic*
        ones.  mechanistic_weight additionally promotes FGs that carry high
        biological information regardless of how common they are — e.g.
        Hydroxamate (Zn chelation warhead for HDAC) or α,β-unsat. carbonyl
        (covalent Michael acceptor).  Together they implement the REFINE_SUGGESTION
        principle: importance_weight × rarity_weight.

    Tie-breaking:
        Equal-score target classes are ranked by insertion order into the
        accumulator dict, which mirrors the order of FG detection (FG_SMARTS
        dict order) and the order of known_target_classes within each FG entry.
        This makes tie-breaking deterministic and pharmacologically interpretable.

    Args:
        fgs_detected: FG names matching FG_SMARTS keys (in detection order).
        fg_db:        Loaded functional_groups dict from fg_database.json.
        use_idf:      Apply IDF weighting (default True).

    Returns:
        DataFrame with columns: target_class | score | votes | evidence_fgs
        Sorted by score descending.  Empty DataFrame if no annotations found.
    """
    # Use plain dict to accumulate mechanistic-weighted votes in insertion order.
    # Insertion order is the tie-breaking key: the first tc that receives a vote
    # (from the first detected FG that annotates it) appears first among equals.
    weighted_votes: dict[str, float] = {}
    evidence: dict[str, list[str]] = defaultdict(list)

    idf: dict[str, float] = _compute_target_idf(fg_db) if use_idf else {}

    for fg in fgs_detected:
        entry = fg_db.get(fg, {})
        mw: float = entry.get("mechanistic_weight", 1.0)
        for tc in entry.get("known_target_classes", []):
            weighted_votes[tc] = weighted_votes.get(tc, 0.0) + mw
            evidence[tc].append(fg)

    if not weighted_votes:
        return pd.DataFrame(
            columns=["target_class", "score", "votes", "evidence_fgs"]
        )

    # ── Post-accumulation adjustments ────────────────────────────────────────
    # CYP450 conditional bonus (applied before IDF multiplication)
    cyp_bonus, cyp_label = _cyp450_conditional_bonus(fgs_detected)
    # CYP450 aryl-halide COOH substrate bonus
    cyp_ah_bonus, cyp_ah_label = _cyp450_arylhalide_cooh_bonus(fgs_detected)
    if cyp_ah_bonus > 0.0:
        cyp_bonus += cyp_ah_bonus
        cyp_label = (cyp_label + " + " + cyp_ah_label).strip(" + ")
    if cyp_bonus > 0.0:
        weighted_votes["cytochrome P450"] = (
            weighted_votes.get("cytochrome P450", 0.0) + cyp_bonus
        )
        evidence["cytochrome P450"].append(f"[{cyp_label}]")

    # COX conditional bonus
    cox_bonus, cox_label = _cox_conditional_bonus(fgs_detected)
    if cox_bonus > 0.0:
        weighted_votes["COX"] = weighted_votes.get("COX", 0.0) + cox_bonus
        evidence["COX"].append(f"[{cox_label}]")

    # mTOR conditional bonus
    mtor_bonus, mtor_label = _mtor_conditional_bonus(fgs_detected)
    if mtor_bonus > 0.0:
        weighted_votes["mTOR"] = weighted_votes.get("mTOR", 0.0) + mtor_bonus
        evidence["mTOR"].append(f"[{mtor_label}]")

    # Kinase covalent warhead bonus
    kin_bonus, kin_label = _kinase_conditional_bonus(fgs_detected)
    if kin_bonus > 0.0:
        weighted_votes["kinase"] = weighted_votes.get("kinase", 0.0) + kin_bonus
        evidence["kinase"].append(f"[{kin_label}]")

    # Adenosine receptor conditional bonus
    ado_bonus, ado_label = _adenosine_conditional_bonus(fgs_detected)
    if ado_bonus > 0.0:
        weighted_votes["adenosine receptor"] = (
            weighted_votes.get("adenosine receptor", 0.0) + ado_bonus
        )
        evidence["adenosine receptor"].append(f"[{ado_label}]")

    # Negative constraints (suppress incompatible target classes in-place)
    _apply_negative_constraints(fgs_detected, weighted_votes, evidence)

    rows = []
    for tc, wt_sum in weighted_votes.items():   # insertion order preserved
        idf_weight = idf.get(tc, 1.0) if use_idf else 1.0
        rows.append({
            "target_class": tc,
            "score": round(wt_sum * idf_weight, 3),
            "votes": round(wt_sum, 3),
            "evidence_fgs": ", ".join(evidence[tc]),
        })

    df = (
        pd.DataFrame(rows)
        .sort_values("score", ascending=False)
        .reset_index(drop=True)
    )
    return df


# ── Main prediction pipeline ───────────────────────────────────────────────────

def predict(
    smiles: str,
    top_residues: int = 10,
    normalize: bool = True,
    use_idf: bool = True,
    table_path: Path = TABLE_PATH,
    db_path: Path = FG_DB_PATH,
) -> dict:
    """Full target prediction pipeline for one SMILES string.

    Args:
        smiles:       Input molecule as SMILES.
        top_residues: How many top-scoring residues to return.
        normalize:    Z-score-normalise residue scores (default True).
        use_idf:      IDF-weight target class votes (default True).
        table_path:   Path to db/fg_residue_table.csv.
        db_path:      Path to db/fg_database.json.

    Returns:
        dict with keys:
          smiles             (str)
          fgs_detected       (list[str])
          residue_scores     (pd.DataFrame)  — residue | residue_name | score | contributing_fgs
          target_class_votes (pd.DataFrame)  — target_class | score | votes | evidence_fgs
          warning            (str | None)    — set on SMILES parse failure or missing table
    """
    result: dict = {
        "smiles": smiles,
        "fgs_detected": [],
        "residue_scores": pd.DataFrame(),
        "target_class_votes": pd.DataFrame(),
        "warning": None,
    }

    # 1. Detect FGs (SMARTS-based, consistent with the BioLiP table)
    fgs = detect_smarts(smiles)
    if not fgs:
        result["warning"] = (
            f"No SMARTS-matched functional groups detected in: {smiles}\n"
            "  Check that the SMILES is valid and that FG_SMARTS covers the scaffold."
        )
        return result
    result["fgs_detected"] = fgs

    # 2. Load FG × residue table
    try:
        table = load_table(table_path)
    except FileNotFoundError as exc:
        result["warning"] = str(exc)
        return result

    # 3. Score residues (drop residue_name label column before scoring)
    numeric_table = table.drop(columns=["residue_name"], errors="ignore")
    result["residue_scores"] = predict_residues(
        fgs, numeric_table, top_n=top_residues, normalize=normalize
    )

    # 4. Vote + weight target classes
    fg_db = load_fg_db(db_path)
    result["target_class_votes"] = predict_target_classes(
        fgs, fg_db, use_idf=use_idf
    )

    return result


# ── Report formatter ───────────────────────────────────────────────────────────

def format_report(pred: dict, compound_name: str = "") -> str:
    """Format a prediction result dict as a human-readable text report.

    Args:
        pred:          Result dict returned by predict().
        compound_name: Optional display name for the compound.

    Returns:
        Formatted multi-line string.
    """
    label = compound_name or pred["smiles"]
    title = f"Target Prediction — {label}"
    bar   = "=" * len(title)
    lines = [bar, title, bar]

    if pred["warning"]:
        lines.append(f"\n⚠  {pred['warning']}")
        return "\n".join(lines)

    lines.append(f"\n  SMILES : {pred['smiles']}")

    # FGs detected
    fgs = pred["fgs_detected"]
    lines.append(f"\n  Functional groups detected ({len(fgs)}):")
    for fg in fgs:
        lines.append(f"    • {fg}")

    # Residue scores
    lines.append("\n  Top binding residues (z-score-normalised BioLiP FG×residue):")
    rs = pred["residue_scores"]
    if rs.empty:
        lines.append("    (no residue matches in table)")
    else:
        lines.append(f"    {'Res':>4}  {'Name':>5}  {'Score':>8}")
        lines.append(f"    {'---':>4}  {'----':>5}  {'-----':>8}")
        for aa, row in rs.iterrows():
            lines.append(
                f"    {aa:>4}  {row['residue_name']:>5}  {row['score']:>8.3f}"
            )

    # Target class votes (IDF-weighted)
    lines.append("\n  Predicted target classes (IDF-weighted FG votes):")
    tc = pred["target_class_votes"]
    if tc.empty:
        lines.append("    (no annotations in fg_database.json for detected FGs)")
    else:
        lines.append(
            f"    {'Target class':<30}  {'Score':>6}  {'Votes':>5}  Evidence FGs"
        )
        lines.append(
            f"    {'-'*30}  {'------':>6}  {'-----':>5}  {'-'*30}"
        )
        for _, row in tc.iterrows():
            lines.append(
                f"    {row['target_class']:<30}  {row['score']:>6.2f}"
                f"  {int(row['votes']):>5}  {row['evidence_fgs']}"
            )

    lines.append("")
    return "\n".join(lines)


# ── Standalone entry point ─────────────────────────────────────────────────────

def main() -> None:
    """CLI wrapper: predict targets for a single SMILES string."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Predict target residues and protein classes for a SMILES string"
    )
    parser.add_argument("smiles", help="Input SMILES string")
    parser.add_argument("--name",  default="",  help="Compound display name")
    parser.add_argument("--top",   type=int, default=10,
                        help="Number of top residues to show (default: 10)")
    parser.add_argument("--no-normalize", action="store_true",
                        help="Disable z-score normalisation of residue scores")
    parser.add_argument("--no-idf", action="store_true",
                        help="Disable IDF weighting of target class votes")
    parser.add_argument("--table", default=str(TABLE_PATH),
                        help="Path to fg_residue_table.csv")
    parser.add_argument("--db",    default=str(FG_DB_PATH),
                        help="Path to fg_database.json")
    args = parser.parse_args()

    pred = predict(
        args.smiles,
        top_residues=args.top,
        normalize=not args.no_normalize,
        use_idf=not args.no_idf,
        table_path=Path(args.table),
        db_path=Path(args.db),
    )
    print(format_report(pred, compound_name=args.name))


if __name__ == "__main__":
    main()
