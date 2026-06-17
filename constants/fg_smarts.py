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
FG_SMARTS defines 42 patterns.  Including Steroid (Python-detected), the
full set used in prediction is 43 functional groups.

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
    # Morpholine: saturated 6-membered ring with O and N at the 1,4 positions.
    #   The morpholine oxygen is the signature hinge-binding group of
    #   ATP-competitive PI3K / mTOR kinase inhibitors (AZD8055, vistusertib,
    #   dactolisib/BEZ235, gedatolisib), H-bonding to the hinge residue
    #   (Val2240 in mTOR / Val851 in PI3Kα). Scaffold marker only — mTOR voting
    #   is handled by the morpholino-diazine conditional rule in
    #   target_predictor.py (Morpholine alone is promiscuous: gefitinib, etc.).
    #   Hierarchical overlap with Ether (ring O) and Tertiary amine (N-alkyl ring
    #   N) is intentional. Aliphatic SMARTS — does not match aromatic oxazines.
    "Morpholine":               "O1CCNCC1",

    # ── Amines ────────────────────────────────────────────────────────────────
    # Each substitution level has a different pKa, H-bond profile, and
    # binding geometry. Amide / sulfonamide / aromatic N excluded — those
    # are captured by Amide, Sulfonamide, and Imidazole respectively.
    "Primary amine":            "[NX3H2;!$(NC=O);!$(NS=O);!$(Nc)]",
    "Secondary amine":          "[NX3H1;!$(NC=O);!$(NS=O);!$(Nc)]",
    "Tertiary amine":           "[NX3H0;!$(N=*);!$(NC=O);!$(NS=O);!$(Nc)]",
    # Imidazole: mimics histidine; coordinates metal ions. Key in CYP450
    # inhibitors (azole antifungals) and histamine receptor ligands.
    # SMARTS covers BOTH free imidazole (nH) AND N-substituted imidazole
    # (n without H) to match azole antifungals such as ketoconazole,
    # clotrimazole, miconazole — where position-1 N is alkylated.
    "Imidazole":                "c1cnc[nH,n]1",
    # Triazole: 5-membered aromatic ring with 3 N atoms (1,2,4-triazole pattern
    #   n-c-n-c-n).  Covers BOTH N-H (1H-triazole) AND N-substituted triazoles.
    #   Key scaffold in triazole antifungals (fluconazole, voriconazole, itraconazole,
    #   posaconazole) that coordinate the heme iron of fungal CYP51 (lanosterol
    #   14α-demethylase) via the free N of the triazole ring.  Distinct from Imidazole
    #   (2N in ring) and Purine (fused bicyclic with 4N).
    #   SMARTS `n1cncn1` matches the 1,2,4-triazole pattern (n-c-n-c-n connectivity)
    #   and was validated on fluconazole, voriconazole (MATCH) and clotrimazole,
    #   metronidazole, omeprazole, purine compounds (NO MATCH).
    "Triazole":                 "n1cncn1",
    # Thiazole: 5-membered aromatic ring with S (position 1) and N (position 3).
    #   1,3-Thiazole isomer; SMARTS covers BOTH free thiazole (N-unsubstituted)
    #   AND N-substituted variants.  Key scaffold in ritonavir-class CYP3A4
    #   inhibitors (HIV protease inhibitors used as pharmacokinetic boosters) where
    #   the thiazole N coordinates heme Fe(III) analogously to imidazole/triazole.
    #   Also present in meloxicam (COX substrate), sulfathiazole (antibacterial),
    #   and various kinase inhibitors.  Substructure of benzothiazole (also matched).
    #   Validated on ritonavir SMILES (2 matches for the 2 thiazole rings) and
    #   confirmed NO match for imidazole, triazole, thiophene, or oxazole.
    "Thiazole":                 "c1cncs1",
    # Benzimidazole: benzo[d]imidazole — benzene fused to imidazole.
    #   Covers both 1H-benzimidazole and N-substituted variants via [nH,n].
    #   Key scaffold in proton pump inhibitors (omeprazole, lansoprazole),
    #   anthelmintic tubulin binders (mebendazole, albendazole), and some
    #   kinase inhibitors.  Contains an imidazole substructure (also detected
    #   separately by the Imidazole pattern) — hierarchical overlap is intentional.
    #   Validated: matches omeprazole, mebendazole, N-methyl-benzimidazole;
    #   does NOT match purine, indole, benzoxazole, or thiazole.
    "Benzimidazole":            "c1ccc2[nH,n]cnc2c1",
    # Pyrimidine: 1,3-diazine — 6-membered aromatic ring with N at positions 1,3.
    #   Core heteroaromatic of ATP-competitive PI3K/mTOR inhibitors (paired with
    #   a morpholine hinge-binder) and many kinase scaffolds. Scaffold marker
    #   only (target_classes=[]): pyrimidine alone is far too promiscuous to vote.
    #   mTOR voting is gated on the Morpholine + (Pyrimidine | Triazine) combo in
    #   target_predictor.py. Substructure of Purine (fused) — hierarchical overlap
    #   is intentional and harmless because it casts no target votes.
    "Pyrimidine":               "c1cncnc1",
    # Triazine: 1,3,5-triazine — 6-membered aromatic ring with three symmetric N.
    #   Replaces pyrimidine as the hinge-anchoring heteroarene in the
    #   morpholino-triazine PI3K/mTOR class (gedatolisib, PF-04691502). Scaffold
    #   marker only; voting handled by the morpholino-diazine conditional rule.
    "Triazine":                 "c1ncncn1",
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
    # Guanidine: -N-C(=N)-N- (three N on the central C). Arginine side-chain
    #   mimetic; like Benzamidine it complements the Asp-lined S1 pocket of
    #   serine proteases (thrombin / trypsin / factor Xa), but via the
    #   guanidinium rather than an aryl amidine. Distinct from Benzamidine
    #   (2 N on C) and matches both neutral and protonated guanidinium.
    "Guanidine":                "[NX3][CX3](=[NX2])[NX3]",

    # ── Sulfur ────────────────────────────────────────────────────────────────
    # Thiol: nucleophilic, metal-coordinating (Cys active sites)
    # Sulfonamide: strong H-bond donor/acceptor, charged at physiological pH
    # Methylsulfone: -SO2CH3 on aromatic ring. COX-2 selectivity pharmacophore
    #   (celecoxib, rofecoxib, valdecoxib). Binds the hydrophilic side-pocket
    #   unique to COX-2 (Val523→Ile in COX-1 blocks this pocket).
    "Thiol":                    "[SX2H]",
    "Sulfonamide":              "[SX4](=O)(=O)[NX3]",
    "Methylsulfone":            "[CX4H3][SX4](=O)(=O)c",

    # ── Metal-binding warheads ────────────────────────────────────────────────
    # Hydroxamate: RC(=O)NHOH. Zinc-chelating warhead for HDAC inhibition
    #   (vorinostat/SAHA, belinostat, panobinostat, pracinostat). Also present
    #   in MMP/ADAM inhibitors (batimastat). The NH bridges to the Zn2+ via
    #   the C=O and OH oxygens (bidentate chelation). The SMARTS requires the
    #   NH to be present (H1 on N), distinguishing from N-alkyl hydroxamates and
    #   acetohydroxamic acid (simpler molecules).
    "Hydroxamate":              "[CX3](=O)[NX3H][OX2H]",
    # Acylsulfonamide: RC(=O)-NH-S(=O)(=O)-R'. Distinguishes the cryptophycin /
    #   epothilone-class tubulin-binding macrolide pharmacophore (C=O adjacent to
    #   SO2 via NH) from classical RS(=O)2NHR sulfonamides. The NH is the key
    #   difference: not a free amine but an acylamide N bridging C=O and SO2.
    #   Also appears in some β-lactamase inhibitors and acylsulfonamide prodrugs.
    "Acylsulfonamide":         "[CX3](=O)[NX3H][SX4](=O)(=O)",

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
    # Anthraquinone (9,10-dioxoanthracene): planar tricyclic quinone. Defining
    #   DNA-intercalating pharmacophore of the anthracycline topoisomerase-II
    #   poisons (doxorubicin, daunorubicin, epirubicin, idarubicin) and
    #   mitoxantrone. The flat polycyclic system stacks between base pairs while
    #   the drug stabilises the Topo-II–DNA cleavage complex. Specific enough to
    #   vote (matches only the anthracycline core, not isolated quinones).
    "Anthraquinone":            "O=C1c2ccccc2C(=O)c2ccccc21",

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
