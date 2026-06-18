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
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


def _find_backend(name: str) -> Optional[str]:
    """Locate a docking executable: PATH, then <env>/Library/bin, then <env>/bin.

    conda envs that aren't `activate`-d don't have Library/bin on PATH, so check
    sys.prefix explicitly (where conda-forge smina/gnina land on Windows)."""
    found = shutil.which(name)
    if found:
        return found
    for sub in ("Library/bin", "bin", "Scripts"):
        for ext in (".exe", ""):
            cand = Path(sys.prefix) / sub / f"{name}{ext}"
            if cand.exists():
                return str(cand)
    return None

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
    """Extract sorted ['RESNAMERESNUM.Interaction', ...] keys from a ProLIF Fingerprint.

    Version-tolerant: reads the to_dataframe() MultiIndex columns and keeps the
    protein-residue + interaction levels of every detected (truthy) interaction.

    The crystal CHAIN id is stripped from each residue (e.g. "ASP189.A" -> "ASP189").
    Serine proteases share the conserved chymotrypsin numbering, so the S1/triad
    residues (ASP189, SER195, GLY216/219, HIS57…) align across trypsin/thrombin/FXa
    regardless of which chain letter a given crystal assigns. Without stripping, a
    query docked into one receptor (chain A) could never match a reference solved in
    another chain (e.g. thrombin chain H) — the Jaccard would be spuriously 0.
    """
    df = fp.to_dataframe()
    keys: set[str] = set()
    for col in df.columns:
        col_t = col if isinstance(col, tuple) else (col,)
        if len(col_t) < 2:
            continue
        prot_res = str(col_t[-2]).rsplit(".", 1)[0]   # drop chain id
        interaction = str(col_t[-1])
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
            backend = self._backend()
            if backend is None:
                return None            # no docking backend → nothing to dock
            mol = self._embed_3d(smiles)
            if mol is None:
                return None
            for cls in classes:
                best_sim, best_ref = self._best_match(mol, cls, backend)
                # Strict '>': a sim exactly at the threshold maps to score 0.0
                # (_sim_to_score), a useless proposal that can only mislead a
                # NO-ANSWER (empty-votes) case. Observed boundary artefact:
                # caffeine docks into the trypsin S1 at sim≈0.60.
                if best_sim > self.sim_threshold:
                    return FallbackProposal(
                        cls, self._sim_to_score(best_sim), self.name,
                        {"ifp_tanimoto": round(best_sim, 3), "ref_pdb": best_ref},
                    )
        except Exception:  # noqa: BLE001 — fallback must never break predict()
            return None
        return None

    def _best_match(self, mol, target_class: str, backend: str) -> tuple[float, Optional[str]]:
        """Dock the query into EACH reference's own receptor and return the best
        (Jaccard, ref_pdb) over the class's actives.

        Each reference IFP was computed in its own crystal frame (FXa, thrombin,
        trypsin…); docking the query into a SINGLE shared receptor would put it in
        the wrong frame and a non-trypsin peptidomimetic could never match its own
        reference (observed: rivaroxaban self-matched at only 0.33 when forced into
        trypsin 3PTB). So dock per reference receptor — query and reference then
        share a frame. Distinct receptors are docked once each (cached).
        """
        best_sim, best_ref = 0.0, None
        ifp_cache: dict[tuple[str, str], set[str]] = {}
        for active in self._load_reference().get(target_class, {}).get("actives", []):
            pdb_id = active.get("source_pdb")
            lig_res = active.get("ligand_resname")
            ref_keys = set(active.get("ifp_keys", []))
            if not pdb_id or not lig_res or not ref_keys:
                continue
            cache_key = (pdb_id, lig_res)
            if cache_key not in ifp_cache:
                ifp_cache[cache_key] = self._dock_and_ifp(mol, pdb_id, lig_res, backend)
            query = ifp_cache[cache_key]
            union = len(query | ref_keys)
            if union == 0:
                continue
            sim = len(query & ref_keys) / union
            if sim > best_sim:
                best_sim, best_ref = sim, pdb_id
        return best_sim, best_ref

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
    def _protonate_smiles(self, smiles: str) -> str:
        """Protonate at pH 7.4 with the SAME OpenBabel model the reference builder
        uses (``build_prolif_reference._protonate_ligand`` → AddHydrogens(.,.,7.4)),
        so query and reference IFPs share ionization state — amidine/guanidine →
        cationic (Cationic + HBDonor S1 keys), carboxylic acid → anionic. Without
        this the neutral query loses the very Cationic/HBDonor keys that define the
        serine-protease S1 match. Falls back to the input SMILES on any failure
        (keeps the zero-regression invariant when OpenBabel is unavailable).
        """
        try:
            from utils.build_prolif_reference import _setup_openbabel  # lazy
            _setup_openbabel()
            from openbabel import pybel  # lazy
            m = pybel.readstring("smi", smiles)
            m.OBMol.AddHydrogens(False, True, 7.4)  # (polar_only, correct_for_pH, pH)
            out = m.write("smi").strip().split()
            return out[0] if out else smiles
        except Exception:  # noqa: BLE001 — protonation is best-effort
            return smiles

    def _embed_3d(self, smiles: str):
        """SMILES → physiological-pH 3D conformer (protonate @ pH 7.4, ETKDGv3 + MMFF)."""
        from rdkit import Chem
        from rdkit.Chem import AllChem
        mol = Chem.MolFromSmiles(self._protonate_smiles(smiles))
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

    def _backend(self) -> Optional[str]:
        """Resolve the docking executable (preferred → smina → gnina) or None."""
        return (_find_backend(self.docking_backend)
                or _find_backend("smina") or _find_backend("gnina"))

    def _prep_receptor(self, pdb_id: str, lig_res: str) -> tuple[Optional[Path], Optional[Path]]:
        """Return (PDBFixer-prepped receptor PDB, autobox-ligand PDB) for one PDB.

        The receptor is prepared with the SAME PDBFixer pass used to build the
        reference library (`build_prolif_reference._prep_protein_pdbfixer`), so the
        query-side IFP is computed against an identically-prepared protein and is
        comparable to that reference's IFP.
        """
        pdb_path = _PDB_DIR / f"{pdb_id}.pdb"
        if not pdb_path.exists():
            return None, None
        from utils.build_prolif_reference import _prep_protein_pdbfixer  # lazy
        rec = _prep_protein_pdbfixer(pdb_path)
        if rec is None:
            return None, None
        box = _PDB_DIR / f"{pdb_id}_{lig_res}_box.pdb"
        if not box.exists():
            import MDAnalysis as mda
            mda.Universe(str(pdb_path)).select_atoms(f"resname {lig_res}").write(str(box))
        return rec, box

    def _dock_and_ifp(self, mol, pdb_id: str, lig_res: str, backend: str) -> set[str]:
        """Dock the query into one reference receptor and return its IFP key-set.

        Query and reference share the receptor frame (same source_pdb), so the
        returned key-set is directly Jaccard-comparable to that reference's keys.
        Returns an empty set on any failure (docking / IFP) so the caller skips it.
        """
        import subprocess
        import tempfile
        import MDAnalysis as mda
        import prolif as plf
        from rdkit import Chem
        rec, box = self._prep_receptor(pdb_id, lig_res)
        if rec is None:
            return set()
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
                return set()
            poses = [m for m in Chem.SDMolSupplier(str(out), removeHs=False)
                     if m is not None]
            if not poses:
                return set()
            prot = plf.Molecule.from_mda(mda.Universe(str(rec)).select_atoms("protein"))
            lig = plf.Molecule.from_rdkit(poses[0])
            fp = plf.Fingerprint()
            fp.run_from_iterable([lig], prot, progress=False)
            return set(ifp_keys_from_fingerprint(fp))

    def _sim_to_score(self, sim: float) -> float:
        """Map IFP Tanimoto (>= threshold) to an FG-comparable score."""
        span = max(1e-6, 1.0 - self.sim_threshold)
        return self.score_scale * (sim - self.sim_threshold) / span
