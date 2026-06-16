"""3D shape descriptors from a SMILES (PMI-based shape + size metrics).

Embeds a single 3D conformer (RDKit ETKDGv3 + MMFF) and computes normalized
principal-moment-of-inertia ratios and related shape/size descriptors. These
capture molecular *shape* (rod / disc / sphere) and *extent* — orthogonal to the
2D functional-group signal — e.g. to test whether elongated CYP450 substrates
separate from compact GPCR ligands.

Descriptors
-----------
npr1, npr2   Normalized PMI ratios (I1/I3, I2/I3). Shape-triangle coordinates:
             rod≈(0,1), disc≈(0.5,0.5), sphere≈(1,1).
radius_gyration   Overall size/extent.
asphericity       0 = spherical, 1 = linear.
eccentricity      0 = spherical, →1 = elongated.
spherocity        Spherocity index (1 = perfect sphere).
inertial_shape_factor

All are pure functions; 3D embedding can fail (returns None) — callers must guard.
"""

from __future__ import annotations

from typing import Optional

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors3D, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")


def embed_3d(smiles: str, seed: int = 0xF00D) -> Optional[Chem.Mol]:
    """SMILES → single 3D conformer (ETKDGv3 + MMFF, explicit H). None on failure."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) != 0:
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, params) != 0:
            return None
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:  # noqa: BLE001 — keep the embedded (un-minimised) conformer
        pass
    return mol


def compute_shape(smiles: str, seed: int = 0xF00D) -> Optional[dict]:
    """Return a dict of 3D shape/size descriptors for a SMILES, or None on failure."""
    mol = embed_3d(smiles, seed=seed)
    if mol is None or mol.GetNumConformers() == 0:
        return None
    try:
        return {
            "npr1": round(rdMolDescriptors.CalcNPR1(mol), 4),
            "npr2": round(rdMolDescriptors.CalcNPR2(mol), 4),
            "radius_gyration": round(Descriptors3D.RadiusOfGyration(mol), 4),
            "asphericity": round(Descriptors3D.Asphericity(mol), 4),
            "eccentricity": round(Descriptors3D.Eccentricity(mol), 4),
            "spherocity": round(Descriptors3D.SpherocityIndex(mol), 4),
            "inertial_shape_factor": round(Descriptors3D.InertialShapeFactor(mol), 4),
        }
    except Exception:  # noqa: BLE001
        return None


# ── Shape-triangle classification (interpretive helper) ─────────────────────────

def shape_class(npr1: float, npr2: float) -> str:
    """Coarse rod/disc/sphere label from NPR coordinates (nearest triangle vertex)."""
    import math
    verts = {"rod": (0.0, 1.0), "disc": (0.5, 0.5), "sphere": (1.0, 1.0)}
    return min(verts, key=lambda k: math.hypot(npr1 - verts[k][0], npr2 - verts[k][1]))


if __name__ == "__main__":
    import sys
    smi = sys.argv[1] if len(sys.argv) > 1 else "CC(=O)Oc1ccccc1C(=O)O"
    s = compute_shape(smi)
    print(smi)
    print(s)
    if s:
        print("shape:", shape_class(s["npr1"], s["npr2"]))
