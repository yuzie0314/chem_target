"""SMARTS patterns for functional group detection and binding-site analysis.

Design principles
-----------------
1. **Specific**: each FG represents a distinct pharmacophoric / binding mode.
2. **Pharmacophore-relevant**: grounded in how chemists and structural biologists
   classify binding interactions.
3. **Natural-product aware**: covers complex scaffolds common in bioactive natural
   products: indole, methylenedioxy, endoperoxide, lactone, coumarin, flavonoid,
   steroid, macrolide.
4. **Hierarchical overlap allowed**: scaffold-level patterns (e.g. Coumarin) may
   co-occur with lower-level patterns (Lactone, Phenyl ring). Each level provides
   independent pharmacological information (binding mode, target class) and is
   counted separately. No two patterns describe *identical* substructures.

Special note — Steroid
----------------------
Steroid scaffold (6-6-6-5 fused tetracyclic) cannot be reliably encoded as a
SMARTS pattern in RDKit because the ``rN`` primitive is tied to the SSSR and
junction atoms between the 5- and 6-membered rings are classified as ``r5``
only, never ``r5;r6``. Detection is therefore implemented as a Python function
``_detect_steroid_core`` in utils/fg_detector.py.  The name "Steroid" is listed
in the docstring and in fg_database.json but is *absent* from FG_SMARTS.

Current FG count
----------------
FG_SMARTS defines 33 patterns.  Including Steroid (Python-detected), the
full set used in prediction is 34 functional groups.

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
    # Ester (acyclic): H-bond acceptor only, metabolically labile; !R excludes
    #   ring esters (lactones), which are captured separately by "Lactone".
    # Amide: resonance-stabilised, donor + acceptor (drug backbone workhorse)
    # Lactone (cyclic ester): same connectivity as Ester but constrained to a
    #   ring; #6X3 and #8X2 match both aliphatic and aromatic (coumarin-type)
    #   forms. Superset of Coumarin and Macrolide.
    "Carboxylic acid":          "C(=O)[OH]",
    "Ester":                    "[#6][CX3;!R](=O)[OX2H0;!R][#6]",
    "Amide":                    "[CX3](=O)[NX3;H1,H2]",
    "Lactone":                  "[#6X3;R](=O)[#8X2;R]",

    # ── Other carbonyls ───────────────────────────────────────────────────────
    # Ketone / aldehyde: pure H-bond acceptor, different reactivity profile.
    # [#6] matches both aliphatic C and aromatic c.
    # Ketone excludes acid/ester/amide by requiring both neighbors to be carbon.
    "Ketone":                   "[#6][#6X3H0](=O)[#6]",
    "Aldehyde":                 "[CX3H1](=O)",

    # ── Hydroxyl groups ───────────────────────────────────────────────────────
    # Aliphatic OH (pKa ~16): H-bond donor/acceptor, weaker than phenol
    # Phenol (pKa ~10):        stronger donor, can ionise, common in polyphenols
    # Catechol (o-diOH):       metal-chelating motif; key in quercetin, EGCG
    "Hydroxyl":                 "[OX2H;!$(Oc);!$(OC=O)]",
    "Phenol":                   "[OX2H]c",
    "Catechol":                 "c1cc([OX2H])c([OX2H])cc1",

    # ── Ether & related ring oxygens ─────────────────────────────────────────
    # Ether: C-O-C where O is not a carbonyl O, not OH, and not directly bonded
    #   to a C=O carbon (i.e., not an ester/acid/carbonate oxygen).
    # Methylenedioxy: -OCH2O- bridge across two adjacent aromatic carbons.
    #   Present in piperine, safrole. Mechanism-based CYP450/MAO inhibitor motif.
    #   Overlaps with Ether (both oxygens match ether pattern), which is
    #   acceptable: each level provides independent pharmacological information.
    "Ether":                    "[OX2;!$(O=*);!$([OX2H]);!$(OC=O)]([#6])[#6]",
    "Methylenedioxy":           "c1ccc2c(c1)OCO2",

    # ── Amines ────────────────────────────────────────────────────────────────
    # Each substitution level has a different pKa, H-bond profile, and
    # binding geometry. Amide / sulfonamide / aromatic N excluded — those
    # are captured by Amide, Sulfonamide, and Imidazole respectively.
    "Primary amine":            "[NX3H2;!$(NC=O);!$(NS=O);!$(Nc)]",
    "Secondary amine":          "[NX3H1;!$(NC=O);!$(NS=O);!$(Nc)]",
    "Tertiary amine":           "[NX3H0;!$(N=*);!$(NC=O);!$(NS=O);!$(Nc)]",
    # Imidazole: mimics histidine; coordinates metal ions. Key in CYP450
    # inhibitors (azole antifungals) and histamine receptor ligands.
    "Imidazole":                "c1cnc[nH]1",
    # Indole: bicyclic N-H aromatic heterocycle (pyrrole fused to benzene).
    # Defining scaffold of tryptamine alkaloids (serotonin, melatonin,
    # vincristine, psilocybin). Excluded from Secondary amine by !$(Nc).
    # Partially overlaps Phenyl ring (benzene portion still matches).
    "Indole":                   "c1ccc2[nH]ccc2c1",
    # Purine: fused 5+6 bicyclic N-heterocycle (imidazole + pyrimidine).
    #   4 N atoms; no H requirement — matches N-H, N-alkyl, and pyridine-type N.
    #   Scaffold of adenine, guanine, hypoxanthine, xanthine, caffeine, and
    #   many kinase inhibitors / antiviral nucleoside analogs.
    #   Superset of Xanthine.
    "Purine":                   "c1ncc2ncnc2n1",
    # Xanthine: 2,6-dioxopurine — the purine with C=O at positions 2 and 6.
    #   Subset of Purine. Scaffold of caffeine (1,3,7-trimethylxanthine),
    #   theophylline, theobromine, and xanthine itself. Key pharmacophores:
    #   adenosine receptor antagonism (caffeine), PDE inhibition (theophylline),
    #   xanthine oxidase inhibition (allopurinol-adjacent).
    "Xanthine":                 "O=c1nc(=O)c2ncnc2n1",

    # ── Other nitrogen ────────────────────────────────────────────────────────
    "Nitrile":                  "C#N",
    "Nitro":                    "[$([NX3](=O)=O),$([NX3+](=O)[O-])]",
    # Benzamidine: amidino group (-C(=NH)NH2) attached to an aromatic ring.
    #   Arginine-mimicking pharmacophore; electrostatically complements the
    #   Asp/Glu-lined S1 pocket of serine proteases (thrombin, trypsin, factor Xa).
    #   Also present in some antiparasitic drugs. Not a guanidine (only 2 N on C).
    "Benzamidine":              "[NX3H2][CX3](=[NX2H1])c",

    # ── Sulfur ────────────────────────────────────────────────────────────────
    # Thiol: nucleophilic, metal-coordinating (Cys active sites)
    # Sulfonamide: strong H-bond donor/acceptor, charged at physiological pH
    # Methylsulfone: -SO2CH3 on aromatic ring. COX-2 selectivity pharmacophore
    #   (celecoxib, rofecoxib, valdecoxib). Binds the hydrophilic side-pocket
    #   unique to COX-2 (Val523→Ile in COX-1 blocks this pocket).
    "Thiol":                    "[SX2H]",
    "Sulfonamide":              "[SX4](=O)(=O)[NX3]",
    "Methylsulfone":            "[CX4H3][SX4](=O)(=O)c",

    # ── Aromatic scaffolds ────────────────────────────────────────────────────
    # Phenyl ring: pure carbocyclic benzene — π-stacking, hydrophobic contacts.
    # Coumarin (2H-chromen-2-one): benzopyranone with C=O at C2 position.
    #   Uses aromatic SMARTS (RDKit treats coumarin as fully aromatic).
    #   Subset of Lactone; target class: MAO inhibitors, anticoagulants (warfarin),
    #   CYP450 substrates, serine proteases.
    # Chromone (4H-chromen-4-one): benzopyranone with C=O at C4 position (isomer
    #   of Coumarin). Backbone of flavonoids (flavone, quercetin, kaempferol).
    #   Uses aromatic SMARTS. Target class: COX, kinases, nuclear receptors.
    #   NOT a lactone (C=O is a ring ketone, ring O is an ether).
    "Phenyl ring":              "c1ccccc1",
    "Coumarin":                 "O=c1ccc2ccccc2o1",
    "Chromone":                 "O=c1ccoc2ccccc12",

    # ── Halogen & reactive warheads ───────────────────────────────────────────
    "Halogen":                  "[F,Cl,Br,I]",
    # Epoxide: electrophilic warhead, covalent binder (terpenoids, mitomycin)
    "Epoxide":                  "[OX2r3]",
    # Endoperoxide: O-O bond in a ring (cyclic peroxide). Pharmacophore of
    #   artemisinin-class antimalarials; reactive with heme Fe(II).
    "Endoperoxide":             "[OX2r][OX2r]",
    # α,β-Unsaturated carbonyl: Michael acceptor (curcumin, parthenolide, etc.)
    "α,β-unsat. carbonyl":      "[CX3](=O)C=C",

    # ── Large-ring & macrocyclic scaffolds ────────────────────────────────────
    # Macrolide: lactone in a ring of ≥12 atoms. Achieved by excluding r3–r11.
    #   Subset of Lactone.  Target class: mTOR (rapamycin), calcineurin
    #   (tacrolimus), ribosome (erythromycin), immunophilins.
    "Macrolide": (
        "[CX3;!r3;!r4;!r5;!r6;!r7;!r8;!r9;!r10;!r11;R](=O)"
        "[OX2;!r3;!r4;!r5;!r6;!r7;!r8;!r9;!r10;!r11;R]"
    ),

    # NOTE — Steroid scaffold (6-6-6-5 fused tetracyclic) is intentionally
    # absent here; it is detected via the Python function
    # ``_detect_steroid_core`` in utils/fg_detector.py.  All other code paths
    # (db/fg_database.json, BioLiP table, target_predictor.py) treat "Steroid"
    # as a regular FG name — only the detection mechanism differs.
}
