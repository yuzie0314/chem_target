"""3D-fallback layer — interface stubs for the structure-based second opinion.

The FG layer is high-precision but leaves NO-ANSWER and FALSE-POSITIVE misses
(see CLAUDE.md). A 3D layer can help, but only as a *gated fallback* registered
via ``target_predictor.register_fallback_3d`` so it never disturbs the FG layer's
high-confidence, zero-regression core.

This module defines the **interface contract** that every 3D fallback must honour,
plus a concrete-but-stubbed ``ProLIFFallback`` (interaction-fingerprint matching,
the planned fix for serine-protease NO-ANSWER peptidomimetics). The detection
pipeline is not implemented yet — ``__call__`` returns ``None`` ("no override"),
so registering this fallback keeps predictions byte-identical (zero regression).
The merge contract (``build_override``) IS implemented, so once the pipeline is
filled in the re-ranking works end-to-end.

Heavy dependencies (ProLIF, a docking backend, optionally RDKit 3D embedding) are
imported lazily *inside* methods, so importing this module never forces them on
the FG core.

Hook contract (matches ``target_predictor._finalize_with_fallback``)
--------------------------------------------------------------------
    fallback(smiles: str, fg_result: dict) -> dict | None
        None  → no override; the FG result is kept verbatim.
        dict  → keys merged into the result via ``result.update(...)``. By
                convention a non-None return contains at least:
                  "target_class_votes" : re-ranked pd.DataFrame
                  "fallback_source"    : str  (e.g. "ProLIF")
                  "fallback_detail"    : dict (provenance: matched ref, score…)

Reference IFP library data contract (per supported target class)
----------------------------------------------------------------
JSON at ``db/prolif_reference_ifp.json`` (gitignored — rebuild offline):
    {
      "serine protease": {
        "target_pdb": "1OYT",                 # protein used for docking
        "binding_site_residues": ["ASP189", "SER195", "GLY216", ...],
        "actives": [
          {"chembl_id": "...", "smiles": "...",
           "ifp_hex": "<ProLIF bitvector as hex>", "source_pdb": "..."},
          ...
        ]
      },
      ...
    }
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd

_ROOT = Path(__file__).resolve().parent.parent
_DEFAULT_REF_IFP = _ROOT / "db" / "prolif_reference_ifp.json"


# ── Proposal + merge contract (shared by all 3D fallbacks) ──────────────────────

@dataclass
class FallbackProposal:
    """A single structure-based target proposal from a 3D fallback."""
    target_class: str
    score: float                      # comparable to FG IDF-weighted scores
    source: str                       # "ProLIF" / "shape" / "Gnina"
    evidence: dict = field(default_factory=dict)   # provenance (ref ligand, sim…)


def build_override(fg_result: dict, proposal: FallbackProposal) -> dict:
    """Build the override dict that re-ranks the FG votes with a 3D proposal.

    Inserts (or boosts) ``proposal.target_class`` in the existing
    ``target_class_votes`` DataFrame at ``proposal.score`` and re-sorts, so a
    confident 3D call can outrank a weak/absent FG winner. Returns the dict the
    hook merges into the result. (This is implemented now; the pipeline that
    *produces* a proposal is the stubbed part.)
    """
    votes = fg_result.get("target_class_votes")
    if votes is None or not isinstance(votes, pd.DataFrame):
        votes = pd.DataFrame(columns=["target_class", "score", "votes", "evidence_fgs"])

    votes = votes.copy()
    tag = f"[3D:{proposal.source}]"
    mask = votes["target_class"] == proposal.target_class
    if mask.any():
        # Boost the existing (under-ranked) class to the 3D score if higher.
        idx = votes.index[mask][0]
        if proposal.score > float(votes.at[idx, "score"]):
            votes.at[idx, "score"] = round(proposal.score, 3)
        votes.at[idx, "evidence_fgs"] = (
            f"{votes.at[idx, 'evidence_fgs']}, {tag}".strip(", ")
        )
    else:
        votes = pd.concat([votes, pd.DataFrame([{
            "target_class": proposal.target_class,
            "score": round(proposal.score, 3),
            "votes": 0.0,
            "evidence_fgs": tag,
        }])], ignore_index=True)

    votes = votes.sort_values("score", ascending=False).reset_index(drop=True)
    return {
        "target_class_votes": votes,
        "fallback_source": proposal.source,
        "fallback_detail": {"target_class": proposal.target_class,
                            "score": proposal.score, **proposal.evidence},
    }


# ── Base interface ──────────────────────────────────────────────────────────────

class Fallback3D:
    """Interface every registered 3D fallback must implement.

    Subclasses set ``name`` and ``supported_classes`` and implement
    ``propose()``; ``__call__`` wires proposal → override and enforces the
    gate-scope contract.
    """
    name: str = "base"
    supported_classes: frozenset[str] = frozenset()

    def propose(self, smiles: str, fg_result: dict) -> Optional[FallbackProposal]:
        """Return a FallbackProposal or None. Override in subclasses."""
        raise NotImplementedError

    def __call__(self, smiles: str, fg_result: dict) -> Optional[dict]:
        proposal = self.propose(smiles, fg_result)
        if proposal is None:
            return None
        return build_override(fg_result, proposal)


# ── ProLIF interaction-fingerprint fallback (STUB) ──────────────────────────────

class ProLIFFallback(Fallback3D):
    """Interaction-fingerprint fallback — planned fix for serine-protease misses.

    Pipeline (each step stubbed; fill in to activate):
      1. _embed_3d(smiles)            RDKit ETKDGv3 + MMFF -> 3D conformer.
      2. _dock(mol, target_class)     dock into the class's reference protein
                                      (smina/gnina/vina) -> best pose.
      3. _compute_ifp(pose, cls)      ProLIF Fingerprint vs binding-site residues
                                      -> interaction bitvector.
      4. _match_reference(ifp, cls)   max Tanimoto to the class's reference-active
                                      IFPs -> (best_sim, best_ref).
      5. if best_sim >= sim_threshold -> FallbackProposal(class, score, ...).

    Until the pipeline is implemented, ``propose`` returns None -> no override ->
    zero regression. Construct with a reference-IFP library and register via
    ``target_predictor.register_fallback_3d(ProLIFFallback(...))``.
    """
    name = "ProLIF"

    def __init__(
        self,
        reference_ifp_path: Path = _DEFAULT_REF_IFP,
        supported_classes: frozenset[str] = frozenset({"serine protease"}),
        sim_threshold: float = 0.6,      # IFP Tanimoto cutoff to propose
        score_scale: float = 8.0,        # maps sim∈[thr,1] -> FG-comparable score
        docking_backend: str = "smina",
    ) -> None:
        self.reference_ifp_path = Path(reference_ifp_path)
        self.supported_classes = frozenset(supported_classes)
        self.sim_threshold = sim_threshold
        self.score_scale = score_scale
        self.docking_backend = docking_backend
        self._ref: Optional[dict] = None   # lazy-loaded IFP library

    # ── public proposal entry point ─────────────────────────────────────────
    def propose(self, smiles: str, fg_result: dict) -> Optional[FallbackProposal]:
        """STUB: returns None (no override). Real flow outlined below."""
        # Scope guard: only attempt classes we hold a reference library for.
        # (The confidence gate in target_predictor already restricts to low/none.)
        # ---- real pipeline (to implement) ----
        # mol  = self._embed_3d(smiles)
        # for cls in self.supported_classes:
        #     pose          = self._dock(mol, cls)
        #     ifp           = self._compute_ifp(pose, cls)
        #     sim, ref      = self._match_reference(ifp, cls)
        #     if sim >= self.sim_threshold:
        #         score = self._sim_to_score(sim)
        #         return FallbackProposal(cls, score, self.name,
        #                                 {"ifp_tanimoto": sim, "ref": ref})
        return None

    # ── reference library ────────────────────────────────────────────────────
    def _load_reference(self) -> dict:
        """Lazy-load db/prolif_reference_ifp.json (see module data contract)."""
        if self._ref is None:
            import json
            if not self.reference_ifp_path.exists():
                self._ref = {}
            else:
                self._ref = json.loads(
                    self.reference_ifp_path.read_text(encoding="utf-8")
                )
        return self._ref

    # ── pipeline steps (STUBS — lazy heavy imports go inside) ─────────────────
    def _embed_3d(self, smiles: str):
        """RDKit ETKDGv3 + MMFF 3D conformer. (to implement)"""
        raise NotImplementedError("ProLIF 3D embedding not yet implemented")

    def _dock(self, mol, target_class: str):
        """Dock into the class's reference protein via self.docking_backend."""
        raise NotImplementedError("ProLIF docking backend not yet implemented")

    def _compute_ifp(self, pose, target_class: str):
        """ProLIF Fingerprint of pose vs binding-site residues -> bitvector."""
        raise NotImplementedError("ProLIF IFP computation not yet implemented")

    def _match_reference(self, ifp, target_class: str) -> tuple[float, Optional[str]]:
        """Max IFP-Tanimoto to the class's reference actives -> (sim, ref_id)."""
        raise NotImplementedError("ProLIF reference matching not yet implemented")

    def _sim_to_score(self, sim: float) -> float:
        """Map IFP Tanimoto (>= threshold) to an FG-comparable score."""
        span = max(1e-6, 1.0 - self.sim_threshold)
        return self.score_scale * (sim - self.sim_threshold) / span
