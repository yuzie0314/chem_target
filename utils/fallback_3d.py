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
        "binding_site_residues": ["ASP189", "SER195", "GLY216", ...],
        "docking": {"receptor_pdb": "3PTB", "autobox_ligand_resname": "BEN"},
        "actives": [
          {"source_pdb": "3PTB", "ligand_resname": "BEN", "note": "...",
           "ifp_keys": ["ASP189.HBDonor", "SER195.HBAcceptor", ...]},
          ...
        ]
      },
      ...
    }
"""

from __future__ import annotations

import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd

_ROOT = Path(__file__).resolve().parent.parent
_DEFAULT_REF_IFP = _ROOT / "db" / "prolif_reference_ifp.json"
_PDB_DIR = _ROOT / "db" / "prolif_pdb"   # cached PDBs (shared with build script)


# ── IFP representation (shared by builder and fallback) ─────────────────────────
# An interaction fingerprint is stored as a SET of "PROTRES.Interaction" keys
# (e.g. "ASP189.HBDonor", "SER195.HBAcceptor") rather than a raw bitvector,
# because ProLIF's bit indexing is not aligned across molecules/runs.  Key-sets
# give a stable, human-readable space and an unambiguous Tanimoto (Jaccard).

def ifp_keys_from_fingerprint(fp) -> list[str]:
    """Extract sorted ['PROTRES.Interaction', ...] keys from a run ProLIF Fingerprint.

    Version-tolerant: reads the to_dataframe() MultiIndex columns and keeps the
    protein-residue + interaction levels of every detected (truthy) interaction.
    """
    df = fp.to_dataframe()
    keys: set[str] = set()
    for col in df.columns:
        col_t = col if isinstance(col, tuple) else (col,)
        if len(col_t) < 2:
            continue
        prot_res, interaction = str(col_t[-2]), str(col_t[-1])
        try:
            on = bool(df[col].to_numpy().any())
        except Exception:  # noqa: BLE001
            on = True
        if on:
            keys.add(f"{prot_res}.{interaction}")
    return sorted(keys)


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
        """Embed → dock → IFP → reference-match → propose (or None).

        Any failure (missing deps, docking error, no reference) is swallowed and
        returns None, so the FG result is kept verbatim — the zero-regression
        invariant holds even if this is registered without ProLIF/docking present.
        """
        try:
            classes = [c for c in self.supported_classes if self._has_reference(c)]
            if not classes:
                return None            # no reference library yet → nothing to match
            mol = self._embed_3d(smiles)
            if mol is None:
                return None
            for cls in classes:
                pose = self._dock(mol, cls)
                if pose is None:
                    continue
                ifp = self._compute_ifp(pose, cls)
                if not ifp:
                    continue
                sim, ref = self._match_reference(ifp, cls)
                if sim >= self.sim_threshold:
                    return FallbackProposal(
                        cls, self._sim_to_score(sim), self.name,
                        {"ifp_tanimoto": round(sim, 3), "ref_pdb": ref},
                    )
        except Exception:  # noqa: BLE001 — fallback must never break predict()
            return None
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

    def _has_reference(self, target_class: str) -> bool:
        return bool(self._load_reference().get(target_class, {}).get("actives"))

    # ── pipeline steps (lazy heavy imports inside each) ───────────────────────
    def _embed_3d(self, smiles: str):
        """SMILES → 3D conformer (RDKit ETKDGv3 + MMFF, explicit H)."""
        from rdkit import Chem
        from rdkit.Chem import AllChem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 0xC0FFEE
        if AllChem.EmbedMolecule(mol, params) != 0:
            params.useRandomCoords = True
            if AllChem.EmbedMolecule(mol, params) != 0:
                return None
        try:
            AllChem.MMFFOptimizeMolecule(mol)
        except Exception:  # noqa: BLE001 — keep the embedded (un-minimised) pose
            pass
        return mol

    def _prep_receptor(self, target_class: str) -> tuple[Optional[Path], Optional[Path]]:
        """Split the cached reference PDB into receptor + autobox-ligand files."""
        dock = self._load_reference().get(target_class, {}).get("docking", {})
        pdb_id = dock.get("receptor_pdb")
        lig_res = dock.get("autobox_ligand_resname")
        if not pdb_id or not lig_res:
            return None, None
        pdb_path = _PDB_DIR / f"{pdb_id}.pdb"
        if not pdb_path.exists():
            return None, None
        rec = _PDB_DIR / f"{pdb_id}_receptor.pdb"
        box = _PDB_DIR / f"{pdb_id}_{lig_res}_box.pdb"
        if not (rec.exists() and box.exists()):
            import MDAnalysis as mda
            u = mda.Universe(str(pdb_path))
            u.select_atoms("protein").write(str(rec))
            u.select_atoms(f"resname {lig_res}").write(str(box))
        return rec, box

    def _dock(self, mol, target_class: str):
        """Dock the query into the class's reference receptor → best-pose RDKit mol."""
        import subprocess
        import tempfile
        from rdkit import Chem
        rec, box = self._prep_receptor(target_class)
        if rec is None:
            return None
        backend = (shutil.which(self.docking_backend)
                   or shutil.which("smina") or shutil.which("gnina"))
        if backend is None:
            return None
        with tempfile.TemporaryDirectory() as td:
            lig_in = Path(td) / "lig.sdf"
            out = Path(td) / "out.sdf"
            w = Chem.SDWriter(str(lig_in)); w.write(mol); w.close()
            cmd = [backend, "--receptor", str(rec), "--ligand", str(lig_in),
                   "--autobox_ligand", str(box), "--autobox_add", "4",
                   "--num_modes", "1", "--exhaustiveness", "8",
                   "--seed", "0",          # deterministic docking
                   "-o", str(out)]
            subprocess.run(cmd, check=True, capture_output=True, timeout=900)
            if not out.exists():
                return None
            poses = [m for m in Chem.SDMolSupplier(str(out), removeHs=False)
                     if m is not None]
            return poses[0] if poses else None

    def _compute_ifp(self, pose, target_class: str) -> set[str]:
        """ProLIF interaction fingerprint of the docked pose vs the receptor."""
        import MDAnalysis as mda
        import prolif as plf
        rec, _ = self._prep_receptor(target_class)
        if rec is None:
            return set()
        prot = plf.Molecule.from_mda(mda.Universe(str(rec)).select_atoms("protein"))
        lig = plf.Molecule.from_rdkit(pose)
        fp = plf.Fingerprint()
        fp.run_from_iterable([lig], prot, progress=False)
        return set(ifp_keys_from_fingerprint(fp))

    def _match_reference(self, ifp_keys, target_class: str) -> tuple[float, Optional[str]]:
        """Max IFP Tanimoto (Jaccard) to the class's reference actives → (sim, ref_pdb)."""
        query = set(ifp_keys)
        best_sim, best_ref = 0.0, None
        for active in self._load_reference().get(target_class, {}).get("actives", []):
            ref = set(active.get("ifp_keys", []))
            union = len(query | ref)
            if union == 0:
                continue
            sim = len(query & ref) / union
            if sim > best_sim:
                best_sim, best_ref = sim, active.get("source_pdb")
        return best_sim, best_ref

    def _sim_to_score(self, sim: float) -> float:
        """Map IFP Tanimoto (>= threshold) to an FG-comparable score."""
        span = max(1e-6, 1.0 - self.sim_threshold)
        return self.score_scale * (sim - self.sim_threshold) / span
