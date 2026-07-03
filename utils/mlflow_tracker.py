"""MLflow experiment tracking for chem_target benchmark runs (opt-in).

Pure side-effect helper: the benchmark pipeline computes its metrics exactly as
before; this module only *records* them to a local MLflow store — a SQLite
backend (``<repo>/mlflow.db``) with artifacts under ``<repo>/mlartifacts``
(MLflow 3 retired the bare file store).  Importing this module has no effect on
prediction logic, and ``mlflow`` is imported lazily so the dependency is only
required when tracking is actually requested (``run_benchmark.py --mlflow``).

Two entry points
----------------
log_benchmark_run   live logging from generate_report() for the current HEAD.
log_historical_run  backfill a single past milestone (known counts only) as a
                    run stamped with the original commit date + SHA, tagged
                    ``reconstructed`` so it is never confused with a real rerun.
"""
from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Optional

EXPERIMENT_NAME = "chem_target-benchmark"

# Core mechanistic classes behind the headline "x/220" number.  The five
# blind-spot classes (COMT, MAO, cysteine protease, topoisomerase,
# xanthine oxidase) are excluded — they are pChEMBL-bias-capped extensions.
CORE_11_CLASSES: frozenset[str] = frozenset({
    "GPCR", "HDAC", "carbonic anhydrase", "tubulin", "nuclear receptor",
    "serine protease", "COX", "kinase", "CYP450", "adenosine receptor", "mTOR",
})


def _repo_root() -> Path:
    """Return the git repo root (falls back to this file's grandparent)."""
    try:
        out = subprocess.run(
            ["git", "rev-parse", "--show-toplevel"],
            capture_output=True, text=True, check=True,
        )
        return Path(out.stdout.strip())
    except (subprocess.CalledProcessError, FileNotFoundError):
        return Path(__file__).resolve().parent.parent


def _tracking_uri() -> str:
    """SQLite tracking backend at <repo>/mlflow.db (MLflow 3 file store retired)."""
    return f"sqlite:///{(_repo_root() / 'mlflow.db').as_posix()}"


def _artifact_uri() -> str:
    """Local artifact root at <repo>/mlartifacts (POSIX-style URI)."""
    return (_repo_root() / "mlartifacts").as_uri()


def _ensure_experiment(client, name: str) -> str:
    """Get-or-create the experiment with our artifact root; return its id."""
    exp = client.get_experiment_by_name(name)
    if exp is not None:
        return exp.experiment_id
    return client.create_experiment(name, artifact_location=_artifact_uri())


def _git_sha(short: bool = True) -> str:
    """Current commit SHA, or 'unknown' outside a git checkout."""
    fmt = ["--short"] if short else []
    try:
        out = subprocess.run(
            ["git", "rev-parse", *fmt, "HEAD"],
            capture_output=True, text=True, check=True,
        )
        return out.stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unknown"


def _git_dirty() -> bool:
    """True if the working tree has uncommitted changes."""
    try:
        out = subprocess.run(
            ["git", "status", "--porcelain"],
            capture_output=True, text=True, check=True,
        )
        return bool(out.stdout.strip())
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False


def _core_subset(cls_df) -> dict[str, float]:
    """Aggregate Top-1/Top-3 over the core-11 mechanistic classes only."""
    core = cls_df[cls_df["target_class"].isin(CORE_11_CLASSES)]
    n = int(core["n_compounds"].sum())
    top1 = int(core["top1_count"].sum())
    top3 = int(core["top3_count"].sum())
    return {
        "core11_n":         n,
        "core11_top1_count": top1,
        "core11_top3_count": top3,
        "core11_top1_acc":  round(top1 / n, 4) if n else 0.0,
        "core11_top3_acc":  round(top3 / n, 4) if n else 0.0,
    }


def log_benchmark_run(
    *,
    mode: str,
    n_all: int,
    n_valid_fg: int,
    top1_all: int,
    top3_all: int,
    mrr_all: float,
    macro_f1: float,
    weighted_f1: float,
    cls_df,
    report_path: Optional[Path] = None,
    summary_path: Optional[Path] = None,
    extra_params: Optional[dict] = None,
) -> Optional[str]:
    """Log one live benchmark run to MLflow; return the run_id (or None on failure).

    All values are the already-computed outputs of generate_report(); this
    function performs no scoring and never raises into the caller.
    """
    try:
        import mlflow
        from mlflow.tracking import MlflowClient
    except ImportError:
        print("[mlflow] mlflow not installed — skipping experiment logging.")
        return None

    mlflow.set_tracking_uri(_tracking_uri())
    _ensure_experiment(MlflowClient(), EXPERIMENT_NAME)
    mlflow.set_experiment(EXPERIMENT_NAME)

    core = _core_subset(cls_df)
    sha = _git_sha()
    dirty = _git_dirty()

    with mlflow.start_run(run_name=f"{mode}-{sha}") as run:
        params = {
            "mode":         mode,
            "git_sha":      sha,
            "git_dirty":    dirty,
            "n_compounds":  n_all,
            "n_classes":    int(cls_df.shape[0]),
        }
        if extra_params:
            params.update(extra_params)
        mlflow.log_params(params)

        mlflow.set_tags({
            "git_sha":       sha,
            "reconstructed": "false",
            "dataset":       f"{n_all}cpd-{cls_df.shape[0]}class",
        })

        mlflow.log_metrics({
            "top1_acc":      round(top1_all / n_all, 4) if n_all else 0.0,
            "top3_acc":      round(top3_all / n_all, 4) if n_all else 0.0,
            "top1_count":    int(top1_all),
            "top3_count":    int(top3_all),
            "n_compounds":   int(n_all),
            "n_valid_fg":    int(n_valid_fg),
            "mrr":           round(float(mrr_all), 4),
            "macro_f1":      round(float(macro_f1), 4),
            "weighted_f1":   round(float(weighted_f1), 4),
            **core,
        })

        # Per-class Top-1 accuracy (namespaced so MLflow groups them).
        for _, r in cls_df.iterrows():
            cls = str(r["target_class"]).replace("/", "-")
            mlflow.log_metric(f"top1_acc.{cls}", float(r["top1_acc"]))

        for p in (report_path, summary_path):
            if p and Path(p).exists():
                mlflow.log_artifact(str(p))

        print(f"[mlflow] logged run {run.info.run_id} "
              f"(core-11 Top-1 {core['core11_top1_count']}/{core['core11_n']})")
        return run.info.run_id


def log_historical_run(
    *,
    run_name: str,
    sha: str,
    date_iso: str,
    n_all: int,
    n_classes: int,
    top1_count: int,
    top3_count: Optional[int],
    note: str = "",
) -> Optional[str]:
    """Backfill one pre-DVC milestone as a reconstructed run.

    ``date_iso`` (e.g. '2026-06-01') stamps the run's start time so the MLflow
    timeline matches git history.  Tagged ``reconstructed=true`` because the
    original dataset snapshot is gone — these are reference points, not reruns.
    """
    try:
        import mlflow
        from mlflow.tracking import MlflowClient
    except ImportError:
        print("[mlflow] mlflow not installed — skipping backfill.")
        return None

    import datetime as _dt

    mlflow.set_tracking_uri(_tracking_uri())
    client = MlflowClient()
    exp_id = _ensure_experiment(client, EXPERIMENT_NAME)

    start_ms = int(
        _dt.datetime.fromisoformat(date_iso)
        .replace(tzinfo=_dt.timezone.utc)
        .timestamp() * 1000
    )
    run = client.create_run(exp_id, start_time=start_ms, run_name=run_name)
    rid = run.info.run_id

    client.log_param(rid, "git_sha", sha)
    client.log_param(rid, "n_compounds", n_all)
    client.log_param(rid, "n_classes", n_classes)
    client.set_tag(rid, "reconstructed", "true")
    client.set_tag(rid, "git_sha", sha)
    client.set_tag(rid, "date", date_iso)
    if note:
        client.set_tag(rid, "note", note)

    top1_acc = round(top1_count / n_all, 4) if n_all else 0.0
    client.log_metric(rid, "top1_acc", top1_acc)
    client.log_metric(rid, "top1_count", top1_count)
    client.log_metric(rid, "n_compounds", n_all)
    # Historically the whole curated benchmark WAS the 11-class core, so these
    # points also populate the core-11 series → one coherent trend with the live
    # run, whose top1_* is on the wider 16-class basis.
    client.log_metric(rid, "core11_top1_acc", top1_acc)
    client.log_metric(rid, "core11_top1_count", top1_count)
    client.log_metric(rid, "core11_n", n_all)
    if top3_count is not None:
        client.log_metric(rid, "top3_acc", round(top3_count / n_all, 4) if n_all else 0.0)
        client.log_metric(rid, "top3_count", top3_count)

    client.set_terminated(rid)
    print(f"[mlflow] backfilled '{run_name}' ({date_iso}, {top1_count}/{n_all})")
    return rid
