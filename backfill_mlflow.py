"""One-off backfill of pre-DVC benchmark milestones into MLflow.

The headline Top-1 percentages recorded in git commit messages are all against
the 11-class / 220-compound *core* curated set — the same population as today's
``core11_top1_acc`` metric — so they form one coherent trend curve.

Because the pre-DVC dataset snapshots were never frozen (benchmark data is
gitignored and was overwritten as classes were added), these counts are
*reconstructed* from the commit-message percentages (pct x 220, rounded) and are
tagged ``reconstructed=true``.  They are reference points, not reruns; the first
fully reproducible run is HEAD logged live via ``run_benchmark.py --mlflow``.

Run once:  python backfill_mlflow.py
Idempotency: re-running creates duplicate runs — only run on a fresh store, or
clear the experiment first.
"""
from __future__ import annotations

from utils.mlflow_tracker import log_historical_run

N_CORE = 220
N_CLASSES = 11

# (sha, date, top1_pct, commit subject) — chronological, 11-class/220 core era.
MILESTONES: list[tuple[str, str, float, str]] = [
    ("d5863c4", "2026-05-31", 59.1, "mechanistic weights baseline"),
    ("d0abc16", "2026-06-01", 62.3, "Ketone->HDAC + Steroid mw + Hydroxamate mw=2.5"),
    ("c347f50", "2026-06-01", 64.1, "COX indole-sulfonamide + CYP450 azole motif"),
    ("89dcc74", "2026-06-01", 64.5, "exclude ab-unsat carbonyl from CYP450 azole rule"),
    ("ba42d7c", "2026-06-01", 65.0, "mTOR macrolide conditional motif"),
    ("15b58e3", "2026-06-01", 65.5, "adenosine receptor Purine motif bonus"),
    ("2ffcb0e", "2026-06-01", 67.7, "kinase covalent ab-unsat carbonyl warhead bonus"),
    ("0ea8ad9", "2026-06-01", 68.2, "kinase Sulfonamide+TertAmine linker bonus"),
    ("406043f", "2026-06-01", 69.1, "CYP450 aryl-halide COOH substrate motif"),
    ("912f4b4", "2026-06-01", 70.0, "CYP450 aryl-halide COOH Rule B (Amide+Ether)"),
    ("3e78a85", "2026-06-01", 70.5, "CYP450 Rule C (ether-amine CYP3A4 scaffold)"),
    ("3372153", "2026-06-01", 70.9, "CYP450 Rule D (amide-halide minimal scaffold)"),
    ("c7dec9d", "2026-06-01", 71.4, "refine CYP450 azole Ketone exclusion"),
    ("49cca41", "2026-06-04", 73.2, "add Thiazole + Benzimidazole SMARTS"),
    ("33eff87", "2026-06-15", 80.5, "morpholino-diazine TORKinib rule (mTOR 5%->85%)"),
    ("e32ec9d", "2026-06-15", 85.5, "pyrimidine ATP-pocket router"),
    ("8a629aa", "2026-06-17", 86.4, "Guanidine FG for serine protease (188->190/220)"),
]


def main() -> None:
    """Log each milestone as a reconstructed historical MLflow run."""
    for sha, date, pct, subject in MILESTONES:
        top1 = round(pct / 100.0 * N_CORE)
        log_historical_run(
            run_name=f"hist-{date}-{sha}",
            sha=sha,
            date_iso=date,
            n_all=N_CORE,
            n_classes=N_CLASSES,
            top1_count=top1,
            top3_count=None,
            note=f"approx from commit msg ({pct}%): {subject}",
        )
    print(f"\nBackfilled {len(MILESTONES)} reconstructed milestones "
          f"(core 11-class / {N_CORE}).")


if __name__ == "__main__":
    main()
