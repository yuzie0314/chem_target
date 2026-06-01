"""Classification metrics for multi-class target prediction evaluation.

Each compound has exactly one true class.  Predictions are a ranked list
(top-K) of class names.  We evaluate with one-vs-rest binary metrics:

  For each class C at threshold K:
    TP  — true class = C  AND  C ∈ top-K predictions
    FP  — true class ≠ C  AND  C ∈ top-K predictions
    FN  — true class = C  AND  C ∉ top-K predictions
    TN  — true class ≠ C  AND  C ∉ top-K predictions

Derived metrics
---------------
  Precision   = TP / (TP+FP)           — of all predictions for C, how many correct?
  Recall      = TP / (TP+FN)           — of all true C, how many retrieved? (Sensitivity)
  Specificity = TN / (TN+FP)           — of all non-C, how many correctly excluded?
  F1          = 2·P·R / (P+R)          — harmonic mean of Precision and Recall
  F2          = 5·P·R / (4·P+R)        — β=2 version; weights missed targets 4× more
                                         → recommended for drug discovery (FN costly)

Aggregation
-----------
  Macro-avg   — unweighted mean across classes (equal voice per class)
  Weighted-avg— support-weighted mean (larger classes contribute more)
"""

from __future__ import annotations

from dataclasses import dataclass, field


# ── Per-class container ──────────────────────────────────────────────────────

@dataclass
class ClassMetrics:
    """One-vs-rest confusion matrix + derived metrics for a single class."""

    name:   str
    n_true: int     # support: number of true examples of this class

    tp: int = field(default=0, repr=False)
    fp: int = field(default=0, repr=False)
    fn: int = field(default=0, repr=False)
    tn: int = field(default=0, repr=False)

    # ── Derived properties ────────────────────────────────────────────────────

    @property
    def precision(self) -> float:
        """TP / (TP+FP).  Returns 0 when no predictions were made for this class."""
        d = self.tp + self.fp
        return self.tp / d if d else 0.0

    @property
    def recall(self) -> float:
        """TP / (TP+FN).  Sensitivity / True Positive Rate."""
        d = self.tp + self.fn
        return self.tp / d if d else 0.0

    @property
    def specificity(self) -> float:
        """TN / (TN+FP).  True Negative Rate."""
        d = self.tn + self.fp
        return self.tn / d if d else 0.0

    @property
    def f1(self) -> float:
        """Harmonic mean of Precision and Recall."""
        p, r = self.precision, self.recall
        d = p + r
        return 2 * p * r / d if d else 0.0

    @property
    def f2(self) -> float:
        """Recall-weighted F-score (β=2).  Penalises false negatives more strongly."""
        p, r = self.precision, self.recall
        d = 4 * p + r
        return 5 * p * r / d if d else 0.0

    def as_dict(self) -> dict:
        """Return all values as a plain dict (useful for CSV rows)."""
        return {
            "class":       self.name,
            "n":           self.n_true,
            "tp":          self.tp,
            "fp":          self.fp,
            "fn":          self.fn,
            "tn":          self.tn,
            "precision":   round(self.precision, 4),
            "recall":      round(self.recall, 4),
            "specificity": round(self.specificity, 4),
            "f1":          round(self.f1, 4),
            "f2":          round(self.f2, 4),
        }


# ── Core computation ─────────────────────────────────────────────────────────

def compute_metrics(
    true_labels:     list[str],
    predicted_sets:  list[list[str | None]],
) -> dict[str, ClassMetrics]:
    """Compute one-vs-rest per-class metrics.

    Args:
        true_labels:    Ground-truth class for each compound (one per compound).
        predicted_sets: For each compound, the list of predicted class labels
                        (up to K elements, highest-ranked first).
                        None / empty strings are silently ignored.

    Returns:
        Mapping from canonical class name → ClassMetrics.
    """
    all_classes = sorted(set(true_labels))
    metrics: dict[str, ClassMetrics] = {
        cls: ClassMetrics(
            name=cls,
            n_true=sum(1 for t in true_labels if t == cls),
        )
        for cls in all_classes
    }

    for true_c, preds in zip(true_labels, predicted_sets):
        # Normalise: deduplicate + drop None / empty
        pred_set: set[str] = {p for p in preds if p}

        for cls, m in metrics.items():
            hit = cls in pred_set
            pos = cls == true_c

            if pos and hit:
                m.tp += 1
            elif pos and not hit:
                m.fn += 1
            elif not pos and hit:
                m.fp += 1
            else:
                m.tn += 1

    return metrics


# ── Aggregation ───────────────────────────────────────────────────────────────

def macro_avg(metrics: dict[str, ClassMetrics]) -> dict[str, float]:
    """Simple (unweighted) mean of each metric across all classes."""
    vals = list(metrics.values())
    if not vals:
        return {}
    n = len(vals)
    return {
        "precision":   sum(m.precision   for m in vals) / n,
        "recall":      sum(m.recall      for m in vals) / n,
        "specificity": sum(m.specificity for m in vals) / n,
        "f1":          sum(m.f1          for m in vals) / n,
        "f2":          sum(m.f2          for m in vals) / n,
    }


def weighted_avg(metrics: dict[str, ClassMetrics]) -> dict[str, float]:
    """Support-weighted mean of each metric across all classes."""
    vals  = list(metrics.values())
    total = sum(m.n_true for m in vals)
    if not total:
        return {}
    return {
        "precision":   sum(m.n_true * m.precision   for m in vals) / total,
        "recall":      sum(m.n_true * m.recall      for m in vals) / total,
        "specificity": sum(m.n_true * m.specificity for m in vals) / total,
        "f1":          sum(m.n_true * m.f1          for m in vals) / total,
        "f2":          sum(m.n_true * m.f2          for m in vals) / total,
    }


# ── Formatted report section ──────────────────────────────────────────────────

def format_metrics_table(
    metrics: dict[str, ClassMetrics],
    title: str = "CLASSIFICATION METRICS (one-vs-rest)",
    show_legend: bool = True,
) -> list[str]:
    """Return a list of formatted plain-text lines for a report section.

    Columns: Class | N | Precision | Recall | Specificity | F1 | F2
    Rows are sorted alphabetically.  Macro and weighted averages appended.
    """
    macro    = macro_avg(metrics)
    weighted = weighted_avg(metrics)

    col_head = (
        f"  {'Class':<24} {'N':>4}  "
        f"{'Prec':>6} {'Recall':>7} {'Spec':>7} {'F1':>6} {'F2':>6}"
    )
    sep = "  " + "-" * 67

    lines = [
        title,
        "-" * 40,
        col_head,
        sep,
    ]

    for cls, m in sorted(metrics.items()):
        lines.append(
            f"  {cls:<24} {m.n_true:>4}  "
            f"{m.precision:>6.3f} {m.recall:>7.3f} {m.specificity:>7.3f} "
            f"{m.f1:>6.3f} {m.f2:>6.3f}"
        )

    lines += [
        sep,
        f"  {'Macro avg':<24} {'':>4}  "
        f"{macro['precision']:>6.3f} {macro['recall']:>7.3f} "
        f"{macro['specificity']:>7.3f} {macro['f1']:>6.3f} {macro['f2']:>6.3f}",
        f"  {'Weighted avg':<24} {'':>4}  "
        f"{weighted['precision']:>6.3f} {weighted['recall']:>7.3f} "
        f"{weighted['specificity']:>7.3f} {weighted['f1']:>6.3f} {weighted['f2']:>6.3f}",
    ]

    if show_legend:
        lines += [
            "",
            "  Definitions (one-vs-rest, binary per class):",
            "    Precision   = TP/(TP+FP)  — when we predict this class, how often correct?",
            "    Recall      = TP/(TP+FN)  — of all true examples, how many retrieved?",
            "    Specificity = TN/(TN+FP)  — of all negatives, how many correctly excluded?",
            "    F1          = 2·P·R/(P+R) — harmonic mean",
            "    F2          = 5·P·R/(4·P+R) — β=2; weighs recall 2× over precision",
            "                  → preferred for drug-discovery (missing a target is costly)",
        ]

    return lines


def format_comparison_metrics_table(
    ct_metrics:  dict[str, ClassMetrics],
    stp_metrics: dict[str, ClassMetrics],
) -> list[str]:
    """Side-by-side F1 / Recall / Precision for chem_target vs STP.

    Intended for the comparison report (short version, no legend).
    """
    all_classes = sorted(set(ct_metrics) | set(stp_metrics))
    ct_macro  = macro_avg(ct_metrics)
    stp_macro = macro_avg(stp_metrics)
    ct_w      = weighted_avg(ct_metrics)
    stp_w     = weighted_avg(stp_metrics)

    col_head = (
        f"  {'Class':<22}  {'N':>4}  "
        f"{'cT-P':>5} {'cT-R':>5} {'cT-F1':>6}  |  "
        f"{'STP-P':>5} {'STP-R':>5} {'STP-F1':>6}"
    )
    sep = "  " + "-" * 76

    lines = [
        "PRECISION / RECALL / F1 COMPARISON (Top-1, one-vs-rest)",
        "-" * 40,
        col_head,
        sep,
    ]

    _z: ClassMetrics = ClassMetrics("", 0)
    for cls in all_classes:
        ct  = ct_metrics.get(cls, _z)
        stp = stp_metrics.get(cls, _z)
        n   = ct.n_true or stp.n_true
        lines.append(
            f"  {cls:<22}  {n:>4}  "
            f"{ct.precision:>5.3f} {ct.recall:>5.3f} {ct.f1:>6.3f}  |  "
            f"{stp.precision:>5.3f} {stp.recall:>5.3f} {stp.f1:>6.3f}"
        )

    lines += [
        sep,
        f"  {'Macro avg':<22}  {'':>4}  "
        f"{ct_macro['precision']:>5.3f} {ct_macro['recall']:>5.3f} "
        f"{ct_macro['f1']:>6.3f}  |  "
        f"{stp_macro['precision']:>5.3f} {stp_macro['recall']:>5.3f} "
        f"{stp_macro['f1']:>6.3f}",
        f"  {'Weighted avg':<22}  {'':>4}  "
        f"{ct_w['precision']:>5.3f} {ct_w['recall']:>5.3f} "
        f"{ct_w['f1']:>6.3f}  |  "
        f"{stp_w['precision']:>5.3f} {stp_w['recall']:>5.3f} "
        f"{stp_w['f1']:>6.3f}",
        "",
        "  P=Precision  R=Recall  F1=harmonic mean(P,R)  [Top-1 predictions only]",
    ]

    return lines
