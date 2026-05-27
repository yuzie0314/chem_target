"""SMARTS patterns for functional group detection and binding-site analysis.

Design principles
-----------------
1. **Specific**: minimal overlap between categories — each FG represents a
   distinct pharmacophoric / binding mode.
2. **Pharmacophore-relevant**: grounded in how chemists and structural biologists
   classify binding interactions.
3. **Natural-product aware**: includes catechol, epoxide, α,β-unsat. carbonyl
   that are common in bioactive natural products.
4. **Non-redundant**: no two patterns describe the same substructure.

Used by
-------
- utils/fg_detector.detect_smarts()    FG profile of a query compound
- utils/interaction_analyzer.py        BioLiP FG × residue co-occurrence table
- utils/target_predictor.py            residue scoring from FG profile
- utils/visualizer.py                  SVG FG highlighting

Rebuilding the BioLiP table
----------------------------
Any change here invalidates db/fg_residue_table.csv.
Rebuild with:
    python utils/interaction_analyzer.py --local db/BioLiP_nr.txt.gz
"""

FG_SMARTS: dict[str, str] = {

    # ── Carboxylic acid & acyl derivatives ────────────────────────────────────
    # COOH: H-bond donor + acceptor + ionisable (pKa ~4.5)
    # Ester: H-bond acceptor only, metabolically labile
    # Amide: resonance-stabilised, donor + acceptor (drug backbone workhorse)
    "Carboxylic acid":          "C(=O)[OH]",
    "Ester":                    "[#6][CX3](=O)[OX2H0][#6]",
    "Amide":                    "[CX3](=O)[NX3;H1,H2]",

    # ── Other carbonyls ───────────────────────────────────────────────────────
    # Ketone / aldehyde: pure H-bond acceptor, different reactivity profile.
    # [#6] matches both aliphatic C and aromatic c, so this correctly captures
    # flavone-type C=O (written as aromatic in some SMILES generators).
    # Excludes acid/ester/amide by requiring both neighbors to be carbon.
    "Ketone":                   "[#6][#6X3H0](=O)[#6]",
    "Aldehyde":                 "[CX3H1](=O)",

    # ── Hydroxyl groups ───────────────────────────────────────────────────────
    # Aliphatic OH (pKa ~16): H-bond donor/acceptor, weaker than phenol
    # Phenol (pKa ~10):        stronger donor, can ionise, common in polyphenols
    # Catechol (o-diOH):       metal-chelating motif; key in quercetin, EGCG
    "Hydroxyl":                 "[OX2H;!$(Oc);!$(OC=O)]",
    "Phenol":                   "[OX2H]c",
    "Catechol":                 "c1cc([OX2H])c([OX2H])cc1",

    # ── Ether ─────────────────────────────────────────────────────────────────
    # C-O-C where O is not a carbonyl O, not OH, and not directly bonded to
    # a C=O carbon (i.e., not an ester/acid/carbonate oxygen).
    "Ether":                    "[OX2;!$(O=*);!$([OX2H]);!$(OC=O)]([#6])[#6]",

    # ── Amines ────────────────────────────────────────────────────────────────
    # Each substitution level has a different pKa, H-bond profile, and
    # binding geometry. Amide / sulfonamide / aromatic N excluded — those
    # are captured by Amide, Sulfonamide, and Imidazole respectively.
    "Primary amine":            "[NX3H2;!$(NC=O);!$(NS=O);!$(Nc)]",
    "Secondary amine":          "[NX3H1;!$(NC=O);!$(NS=O);!$(Nc)]",
    "Tertiary amine":           "[NX3H0;!$(N=*);!$(NC=O);!$(NS=O);!$(Nc)]",
    "Imidazole":                "c1cnc[nH]1",

    # ── Other nitrogen ────────────────────────────────────────────────────────
    "Nitrile":                  "C#N",
    "Nitro":                    "[$([NX3](=O)=O),$([NX3+](=O)[O-])]",

    # ── Sulfur ────────────────────────────────────────────────────────────────
    # Thiol: nucleophilic, metal-coordinating (Cys active sites)
    # Sulfonamide: strong H-bond donor/acceptor, charged at physiological pH
    "Thiol":                    "[SX2H]",
    "Sulfonamide":              "[SX4](=O)(=O)[NX3]",

    # ── Aromatic ──────────────────────────────────────────────────────────────
    # Pure carbocyclic benzene ring — π-stacking, hydrophobic interactions.
    # (Removed duplicate "Benzene" entry that used the same SMARTS.)
    "Phenyl ring":              "c1ccccc1",

    # ── Halogen & reactive warheads ───────────────────────────────────────────
    "Halogen":                  "[F,Cl,Br,I]",
    # Epoxide: electrophilic warhead, covalent binder (terpenoids, mitomycin)
    "Epoxide":                  "[OX2r3]",
    # α,β-Unsaturated carbonyl: Michael acceptor (curcumin, parthenolide, etc.)
    "α,β-unsat. carbonyl":      "[CX3](=O)C=C",
}
