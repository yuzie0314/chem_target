"""Regression safety net for the FG-detection + target-prediction pipeline.

Uses the stdlib ``unittest`` (no extra dependency — the project env has no
pytest) so it runs anywhere the package imports:

    conda run -n chem_target python -m unittest discover -s tests -v

The load-bearing test is ``test_curated_core11_headline``: it locks the headline
190/220 core-11 Top-1 figure against the frozen curated benchmark, turning the
previously-unreproducible manual benchmark into an automated assertion.  Any rule
change that moves the headline will now fail CI-style instead of silently.
"""

from __future__ import annotations

import csv
import unittest
from pathlib import Path

_ROOT = Path(__file__).resolve().parent.parent
import sys  # noqa: E402

sys.path.insert(0, str(_ROOT))

from utils.fg_detector import detect_smarts  # noqa: E402
from utils.target_predictor import (  # noqa: E402
    assess_confidence,
    predict,
    predict_target_classes,
    load_fg_db,
)

CORE11 = {
    "mTOR", "tubulin", "carbonic anhydrase", "HDAC", "CYP450", "COX",
    "serine protease", "adenosine receptor", "GPCR", "kinase", "nuclear receptor",
}
CURATED_CSV = _ROOT / "data" / "benchmark" / "curated" / "compounds.csv"


class TestFGDetection(unittest.TestCase):
    """SMARTS + Python detectors on known molecules."""

    def test_aspirin_functional_groups(self) -> None:
        fgs = set(detect_smarts("CC(=O)Oc1ccccc1C(=O)O"))
        self.assertEqual(fgs, {"Carboxylic acid", "Ester", "Phenyl ring"})

    def test_invalid_smiles_returns_empty(self) -> None:
        self.assertEqual(detect_smarts("not_a_smiles"), [])

    def test_steroid_python_detector(self) -> None:
        # Testosterone — 6-6-6-5 fused tetracyclic scaffold.
        fgs = detect_smarts("CC12CCC3C(CCC4=CC(=O)CCC34C)C1CCC2O")
        self.assertIn("Steroid", fgs)


class TestPredictionInvariants(unittest.TestCase):
    """Pipeline-level invariants that must never silently break."""

    def test_aspirin_prediction_is_high_confidence(self) -> None:
        res = predict("CC(=O)Oc1ccccc1C(=O)O")
        self.assertEqual(res["confidence"], "high")
        # Zero-regression invariant: the default 3D fallback never fires.
        self.assertFalse(res["fallback_applied"])

    def test_no_fg_routes_to_none_confidence(self) -> None:
        res = predict("[He]")  # noble gas: no FG match
        self.assertIn(res["confidence"], {"none", "low"})

    def test_residue_scores_independent_of_target_votes(self) -> None:
        # The BioLiP residue table must not be an input to target-class voting:
        # target votes depend only on fgs_detected × fg_database.
        smiles = "CC(=O)Oc1ccccc1C(=O)O"
        fgs = detect_smarts(smiles)
        votes = predict_target_classes(fgs, load_fg_db())
        res = predict(smiles)
        # Same FG list → identical target-class ranking regardless of residues.
        self.assertEqual(
            list(votes["target_class"]),
            list(res["target_class_votes"]["target_class"]),
        )

    def test_assess_confidence_empty_is_none(self) -> None:
        import pandas as pd
        self.assertEqual(
            assess_confidence(pd.DataFrame(), []), "none",
        )

    def test_aryl_nitrile_not_grabbed_as_cysteine_protease(self) -> None:
        # Enobosarm: an aryl-nitrile + amide androgen (nuclear-receptor) ligand.
        # The cathepsin rule keys on "Non-aryl nitrile", so an inert aryl nitrile
        # must NOT be mis-predicted as cysteine protease.
        res = predict("C[C@](O)(COc1ccc(C#N)cc1)C(=O)Nc1ccc(C#N)c(C(F)(F)F)c1")
        votes = res["target_class_votes"]
        self.assertFalse(votes.empty)
        self.assertNotEqual(str(votes.iloc[0]["target_class"]), "cysteine protease")

    def test_empty_votes_after_negative_constraint_does_not_crash(self) -> None:
        # A fused azolo-diazine + pyrimidine whose only class vote is CYP450 gets
        # that vote suppressed by the negative constraint, emptying the votes.
        # predict_target_classes must return an empty frame, not raise.
        fgs = ["Fused azolo-diazine", "Pyrimidine", "Triazole"]
        df = predict_target_classes(fgs, load_fg_db())
        self.assertIn("target_class", df.columns)


class TestCuratedBenchmarkRegression(unittest.TestCase):
    """Lock the headline accuracy against the frozen tuning set."""

    @classmethod
    def setUpClass(cls) -> None:
        if not CURATED_CSV.exists():
            raise unittest.SkipTest(f"frozen benchmark missing: {CURATED_CSV}")
        with open(CURATED_CSV, encoding="utf-8") as fh:
            cls.rows = [r for r in csv.DictReader(fh)
                        if r["true_target_class"] in CORE11]

    def _canonical(self, predicted: str) -> str:
        # Mirror run_benchmark's alias mapping so the count matches the pipeline.
        from run_benchmark import _map_predicted_class
        return _map_predicted_class(predicted)

    def test_curated_core11_headline(self) -> None:
        hits = 0
        for r in self.rows:
            res = predict(r["smiles"])
            votes = res["target_class_votes"]
            if votes.empty:
                continue
            top1 = self._canonical(str(votes.iloc[0]["target_class"]))
            if top1 == r["true_target_class"]:
                hits += 1
        n = len(self.rows)
        self.assertEqual(
            n, 220,
            f"expected 220 core-11 curated compounds, found {n}",
        )
        # Headline is 185/220 = 84.1% (was 190 before the 2026-07-03 mechanistic
        # de-overfitting: removed the ID-tuned aryl-halide CYP450 rules → −5 CYP450
        # tuning, +0 held-out). Lock it: a rule change that moves this must be
        # deliberate (update the expected value with justification).
        self.assertEqual(
            hits, 185,
            f"core-11 Top-1 regressed/changed: {hits}/220 (expected 185). "
            "If intended, update this assertion and the README/CLAUDE headline.",
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
