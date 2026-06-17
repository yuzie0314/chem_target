\# CLAUDE.md — chem\_target project



\## Who you are working with

Yu-Jen Chang — pharmacist, computational chemist, data scientist.

Deep functional group intuition from pharmacy background.

Primary docking tools: iGEMDOCK, SiMMap (not AutoDock/Glide).

Current stack: Python 3.10, RDKit, OpenBabel, pandas, conda.



\---



\## Project goal

A compound-to-target prediction tool (reverse docking / target fishing).

Input: molecule (SMILES or SDF). Output: predicted enzyme/protein targets

based on functional group analysis. Long-term: consultancy tool + SaaS.



\---



\## Architecture rules — always follow these

```
chem_target/
├── constants/      # Static lookup tables only. No logic.
│   ├── fg_names.py     # FG_NAMES: rdkit fr_* → human readable (legacy)
│   └── fg_smarts.py    # FG_SMARTS: name → SMARTS (42 patterns + Steroid Python = 43 total)
├── db/             # Auto-generated. Never hand-edit.
│   ├── fg_database.json       # 43 FG metadata: smarts/targets/mechanistic_weight
│   ├── fg_residue_table.csv   # 42 SMARTS + Steroid × 20 AA BioLiP co-occurrence matrix
│   ├── ccd_smiles_cache.json  # RCSB CCD SMILES cache
│   └── residue_3d_poses.json  # Cα + ligand centroid 3D records
├── utils/          # Pure functions. No side effects where possible.
│   ├── fg_detector.py         # detect_smarts(), _detect_steroid_core(), _detect_fused_azolo_diazine()
│   ├── target_predictor.py    # IDF × mw scoring + conditional rules + _pyrimidine_router() + confidence gate / register_fallback_3d
│   ├── fallback_3d.py         # 3D-fallback interface (Fallback3D, ProLIFFallback stub, build_override) — lazy heavy deps
│   ├── build_prolif_reference.py  # env check + ProLIF reference-IFP builder (check|build) → db/prolif_reference_ifp.json
│   ├── interaction_analyzer.py  # BioLiP → fg_residue_table.csv
│   ├── pose_extractor.py      # 3D pose extractor → residue_3d_poses.json
│   ├── db_updater.py          # PubChem/ChEMBL → fg_database.json
│   ├── io_handler.py          # CSV/SDF input parsing
│   ├── report_generator.py    # HTML individual + batch reports
│   └── visualizer.py          # RDKit SVG output
├── data/           # User input files
├── output/         # Generated output — CSV, SVG, HTML, reports
├── run_benchmark.py   # 11-class × 20-compound benchmark pipeline
└── main.py         # Entry point only. Thin. No logic here.
```

**Never put constants inside logic files.**
**Never put logic inside main.py.**
**Never hardcode values that belong in constants/ or db/.**
**Never write to db/ manually** — use `interaction_analyzer.py`, `pose_extractor.py`, or `db_updater.py`.



\---



\## Coding conventions



\- Python 3.10

\- Type hints on all function signatures

\- Docstrings on all functions (one-line minimum)

\- Constants in UPPER\_SNAKE\_CASE

\- Functions in snake\_case

\- Early returns over nested if-else

\- Explicit is better than implicit — no magic numbers

\-Record progress using git add/git commit and write down modifications in details and after finish major task, please remember to push the current progress to the remote

\---



\## DB update strategy



\- Primary source: PubChem API (`https://pubchem.ncbi.nlm.nih.gov/rest/pug`)

\- Secondary source: ChEMBL API (`https://www.ebi.ac.uk/chembl/api/data`)

\- DB file: `db/fg\_database.json`

\- Always include `last\_updated` timestamp in DB

\- DB updater must be runnable as standalone: `python utils/db\_updater.py`



\### fg\_database.json schema

```json

{

&#x20; "last\_updated": "YYYY-MM-DD",

&#x20; "sources": \["pubchem", "chembl"],

&#x20; "functional\_groups": {

&#x20;   "Carboxylic acid": {

&#x20;     "smarts": "C(=O)\[OH]",

&#x20;     "pubchem\_cid": 280,

&#x20;     "chembl\_id": null,

&#x20;     "description": "...",

&#x20;     "known\_target\_classes": \["protease", "transporter"]

&#x20;   }

&#x20; }

}

```



\---



\## Input formats (via OpenBabel)



Currently supported (`utils/io_handler.read_file` auto-dispatches by extension):

\- CSV: col1 = compound name, col2 = SMILES (`read_csv`)
\- SDF/.sd: **RDKit `SDMolSupplier`** (native/robust); name from `--name-property` tag → title (`_Name`) → `mol_<i>`; SMILES via RDKit (`read_sdf`)
\- MOL2/.mol: **OpenBabel (pybel)** via `read_mol2`; RDKit's MOL2 reader is weak
\- SMI/.smiles: whitespace-delimited `<SMILES> [name]` per line; header + invalid lines skipped (`read_smiles_file`)
\- InChI/.ich: `InChI=… [name]` per line → RDKit `MolFromInchi` (`read_inchi_file`)

NOTE on OpenBabel: conda's OpenBabel format plugins do NOT load on Windows without a DLL/PATH fix —
`io_handler._setup_openbabel()` (and the same helper in `build_prolif_reference.py`) sets
`BABEL_LIBDIR`/`BABEL_DATADIR` + `os.add_dll_directory` from `sys.prefix`. SDF/SMILES/InChI deliberately
use RDKit (not OpenBabel) to avoid this fragility; OpenBabel is used only for MOL2.

Planned (extend io\_handler.py): InChIKey lookup (needs external resolver — not structural)



Use OpenBabel for all format conversions. Do not write custom parsers.



\---



\## Output formats



\- `fg\_abundance\_table.csv` — rows=functional groups, cols=compounds, values=integer counts

\- `output/images/{compound\_name}\_fg.svg` — molecule with highlighted functional groups

\- Future: PDF/HTML report for non-technical clients



\---



\## Key dependencies



```

rdkit          # Core cheminformatics

openbabel      # Format conversion

pandas         # Tables

requests       # API calls to PubChem / ChEMBL

matplotlib     # Plotting (future)

```



conda environment name: `chem\_target`



\---



\## What NOT to do

\- Do not use AutoDock or Glide — developer uses iGEMDOCK + SiMMap
\- Do not merge constants into logic files
\- Do not write to db/ manually — always via interaction\_analyzer.py / pose\_extractor.py / db\_updater.py
\- Do not over-engineer early — get it working first, then clean up
\- Do not ignore OpenBabel for format handling — it is already installed



\---



\## Current status (2026-06-15)

| Component | Status |
|---|---|
| conda env + rdkit + openbabel | ✅ Done |
| 42 FG SMARTS + Steroid Python = 43 total (`constants/fg_smarts.py`) | ✅ Done |
| `db/fg_database.json` (43 entries incl. …+Pyrimidine+Triazine+Anthraquinone+Guanidine, mechanistic_weight) | ✅ Done |
| Serine protease Guanidine FG (arginine-mimetic → SP, mw 3.0) | ✅ 13/20 (+1, zero regression) |
| `db/fg_residue_table.csv` (BioLiP rebuild, 42 SMARTS + Steroid columns) | ✅ Done |
| `db/residue_3d_poses.json` + `db/local_env/*.sdf` | ✅ Done |
| `utils/target_predictor.py` (IDF × mechanistic_weight) | ✅ Done |
| `utils/report_generator.py` (HTML individual + batch) | ✅ Done |
| `run_benchmark.py` (11-class × 20-compound curated) | ✅ Done |
| **Core 11-class Top-1: 190/220 = 86.4%** (mechanistic classes) | ✅ Current best |
| **Core 11-class Top-3: 197/220 = 89.5%** | ✅ Current best |
| Blind-spot rule-backed: MAO 2/20, COMT 8/20, cysteine protease 12/20, topoisomerase 5/20 | ✅ |
| MAO covalent-warhead rule (Propargylamine/Hydrazine) | ✅ Done |
| COMT (nitrocatechol via existing Phenol+Catechol) | ✅ 8/20 (pChEMBL-bias limited) |
| Cysteine protease nitrile-warhead rule (Nitrile+Amide, gated) | ✅ 12/20 (zero collision) |
| Topoisomerase Anthraquinone voting FG (anthracyclines) | ✅ 5/20 (mw=2.5→topo) |
| CYP450 conditional motif scoring (azole rule, Thiazole added, Pyrimidine guard) | ✅ Done |
| Negative constraint rules (Hydroxamate/Thiol/Acylsulfonamide + fused-azolo-diazine → suppress CYP450) | ✅ Done |
| COX indole-sulfonamide motif | ✅ Done |
| mTOR macrolide conditional motif (rapalog) | ✅ Done |
| **Pyrimidine router** (mutually-exclusive: morpholino→mTOR / fused-azolo→adenosine / mono→kinase) | ✅ Done |
| Adenosine receptor Purine bonus | ✅ Done |
| Kinase α,β-unsat carbonyl covalent warhead bonus | ✅ Done |
| Thiazole SMARTS + BioLiP table rebuild | ✅ Done |
| Benzimidazole SMARTS (scaffold detection only, no target votes) | ✅ Done |
| Morpholine + Pyrimidine + Triazine SMARTS (scaffold markers, no target votes) | ✅ Done |
| `_detect_fused_azolo_diazine` Python detector (routing-only, not in fg_database → no IDF impact) | ✅ Done |
| FG-confidence tiering + 3D-fallback routing skeleton (stub, zero-regression) | ✅ Done (stub) |
| 方案 4 fused-N core: azolo-diazine detector (functional, +7) + Quinazoline/Pyrrolopyrimidine/Pyridopyrimidine/Benzoxazole (annotation-only) | ✅ Done |
| Input formats: CSV + SDF(RDKit) + MOL2(OpenBabel) + SMI + InChI — `read_file` dispatch, wired in main.py fg/predict | ✅ Done |
| 3D shape descriptors (`utils/shape_descriptors.py`: NPR1/NPR2, Rg, asphericity…) — opt-in `predict --shape` | ✅ Done (informational) |

---

\## Benchmark per-class results (current best, branch: dev/validation)

| Class | Top-1 | Top-3 | Notes |
|---|---|---|---|
| GPCR | 20/20 = 100% | 20/20 | ✅ |
| HDAC | 20/20 = 100% | 20/20 | ✅ |
| Carbonic anhydrase | 20/20 = 100% | 20/20 | ✅ |
| Tubulin | 19/20 = 95% | 19/20 | -1 GS-9256 (thiazole+ether → CYP rule; profile indistinguishable from ritonavir-class) |
| Nuclear receptor | 17/20 = 85% | 20/20 | +1 from Guanidine FG IDF shift; remaining: Acylsulfonamide→tubulin + structural |
| Serine protease | 13/20 = 65% | 13/20 | +1 by Guanidine FG (CHEMBL353760, arginine-mimetic→SP). 7 remaining peptidomimetics have NO amidine/guanidine S1 group (structural) |
| COX | 15/20 = 75% | 17/20 | Fixed +4 by Indole+Sulfonamide motif |
| Kinase | 18/20 = 90% | 18/20 | +4 by pyrimidine router (mono-pyrimidine→kinase, branch 3); earlier +6 αβunsat+Sulfonamide. 2 remaining: 1 strong-GPCR (CHEMBL5270693), 1 Steroid |
| CYP450 | 19/20 = 95% | 19/20 | Fixed +12 total; 5 ritonavir-class by Thiazole SMARTS; 1 TAZAROTENIC ACID structural. Pyrimidine guard added (no CYP450 TP has pyrimidine) |
| Adenosine receptor | 12/20 = 60% | 12/20 | +7 by pyrimidine router (fused-azolo-diazine→adenosine, branch 2). 8 remaining: no purine-mimetic core (Phenol/Halogen sparse, or Nitrile/Steroid) |
| mTOR | 17/20 = 85% | 17/20 | Fixed +16 by morpholino-diazine motif (16 ATP-competitive TORKinibs); +SIROLIMUS by macrolide rule. 3 remaining have no morpholine (SAPANISERTIB, CHEMBL3645910, CHEMBL3681183) |
| COMT | 8/20 = 40% | — | nitrocatechol (entacapone/opicapone) via existing Phenol+Catechol; other 12 = research series w/o nitrocatechol (pChEMBL-bias) |
| MAO | 2/20 = 10% | — | propargylamine (clorgiline) via MAO warhead rule; 18 = single research series (Sec/Tert amine, no MAO pharmacophore) |
| cysteine protease | 12/20 = 60% | — | nitrile-warhead cathepsin inhibitors (odanacatib class) via Nitrile+Amide rule; zero collision. 8 remaining lack nitrile |
| topoisomerase | 5/20 = 25% | — | anthracyclines (doxorubicin/daunorubicin/epirubicin/idarubicin/nemorubicin) via Anthraquinone voting FG (mw=2.5, sole topo annotator → IDF≈3.7). 15 = research series w/o intercalator core |
| xanthine oxidase | 0/20 | — | pChEMBL-bias: research series (Phenol+Amide+Pyrimidine), no allopurinol/febuxostat pharmacophore; pyrimidine→kinase collision. Structural limit |

**Note on the 13-class extended set (260 compounds):** MAO + COMT were added 2026-06-15 with
*correct* ChEMBL target IDs (MAO=CHEMBL1951/2039, COMT=CHEMBL2023). Their low scores reflect a
**pChEMBL-sampling bias**: the top-20 highest-affinity ChEMBL compounds for these targets are modern
research analogs that lack the classic covalent pharmacophores (propargylamine, nitrocatechol) the
marketed drugs carry — not a rule gap. Adding the warhead rule cannot capture analogs that don't
carry the warhead. The other 5 blind-spot classes (PDE/topoisomerase/ribosome/XO/cysteine protease)
remain un-added; see Next tasks.

✅ **Data bug FIXED (2026-06-15):** `run_stp_comparison.py` previously had a local `_TARGET_CLASS_MAP`
with *wrong* ChEMBL IDs for the 7 blind-spot classes (e.g. MAO→CHEMBL2828 returned DARUNAVIR,
topoisomerase→CHEMBL3952 returned the opioid JDTIC). It now `from run_benchmark import TARGET_CLASS_MAP`
— **single source of truth**, so the two can never diverge again. `generate_comparison` also restricts
to the true compound intersection (`restrict_names`), keeping the head-to-head fair after benchmark
expansion. NOTE: the cached `stp_raw.csv` still holds the old mislabelled STP rows for those 7 classes
(fetched with the bad map); re-query STP to refresh them — but the published 220/11-class fair
comparison never used them, so it is unaffected.

---

\## mechanistic_weight assignments (fg_database.json)

| FG | mw | Target class | Rationale |
|---|---|---|---|
| Benzamidine | 3.0 | serine protease | S1-pocket bidentate H-bond |
| Guanidine | 3.0 | serine protease | Arginine-mimetic guanidinium → S1 pocket (parallel to Benzamidine) |
| Hydroxamate | 2.5 | HDAC | Bidentate Zn chelation |
| Sulfonamide | 2.0 | carbonic anhydrase | Zn coordination |
| Acylsulfonamide | 2.0 | tubulin | Macrolide warhead |
| Ketone | 2.0 | HDAC | α-keto warhead in HDAC natural products |
| Steroid | 2.0 | nuclear receptor (+ subtypes) | Steroidal scaffold → NR |
| Anthraquinone | 2.5 | topoisomerase | Planar tricyclic quinone DNA-intercalator (anthracyclines/mitoxantrone); sole topo annotator → high IDF |
| Triazole | 1.5 | cytochrome P450 | Triazole antifungal heme-Fe coordination (fluconazole-class) |
| Thiazole | 1.5 | cytochrome P450 | Ritonavir-class CYP3A4 inhibitor; thiazole N coordinates heme Fe analogously to imidazole |
| Benzimidazole | 1.0 | (none — scaffold marker) | Benzene+imidazole fused; Imidazole FG already handles CYP450 voting for benzimidazole compounds |
| Morpholine | 1.0 | (none — scaffold marker) | PI3K/mTOR hinge-binder; promiscuous alone (gefitinib). Voting via pyrimidine router (branch 1) |
| Pyrimidine | 1.0 | (none — scaffold marker) | 1,3-diazine hinge-anchor; promiscuous alone. Voting via pyrimidine router (branches 1–3) |
| Triazine | 1.0 | (none — scaffold marker) | 1,3,5-triazine hinge-anchor (gedatolisib class). Voting via pyrimidine router (branch 1) |
| All others | 1.0 | — | Default |

**Routing-only Python detector (NOT in fg_database.json → no IDF impact):**
`Fused azolo-diazine` — `utils/fg_detector._detect_fused_azolo_diazine`: aromatic 5-ring(≥2N) fused to 6-ring(≥2N). Purine / triazolopyrimidine / pyrazolopyrimidine core. Consumed only by the pyrimidine router (branch 2 → adenosine) and the CYP450 negative constraint. Appears in `fgs_detected` but casts no votes.

---

\## Conditional scoring rules (utils/target_predictor.py)

All rules are pre-IDF bonuses (multiplied by IDF before adding to final score).

| Rule | Condition | Target | Bonus | Rationale |
|---|---|---|---|---|
| CYP450 azole | **free heme azole** (`_has_free_heme_azole`) + {Phenyl/Ether/Halogen}, Ketone only if no Amide/TertAmine, no Purine/αβunsat/Sulfonamide | cytochrome P450 | +2.0 | Azole heme-Fe coordination (fluconazole/voriconazole/ketoconazole/ritonavir class). Free azole = Triazole OR Thiazole OR (Imidazole & not Benzimidazole), AND not Fused-azolo-diazine — excludes purine-mimetic fused cores (adenosine/kinase) while keeping voriconazole (free triazole + separate fluoropyrimidine) |
| CYP450 aryl-COOH A | COOH + Phenyl + Halogen, no Amide, no Ether | cytochrome P450 | +1.5 | Minimal aryl-halide CYP substrate |
| CYP450 aryl-COOH B | COOH + Amide + Ether + Phenyl + Halogen | cytochrome P450 | +1.5 | Extended aryl-halide CYP substrate |
| CYP450 ether-amine | Ether + TertAmine + Phenyl + Halogen, no Lactone/Amide/Nitrile | cytochrome P450 | +1.5 | CYP3A4 scaffold (aprepitant-type) |
| CYP450 amide-halide | Amide + Phenyl + Halogen, no Sulfonamide/COOH/Imidazole/αβunsat/Ether | cytochrome P450 | +0.6 | Minimal amide-halide CYP substrate (raised from 0.5 to compensate IDF shift from Thiazole) |
| COX indole-sulfonamide | Indole + Sulfonamide | COX | +2.0 | Indole scaffold + COX-2 selectivity pocket |
| mTOR macrolide | Macrolide, no Thiol/αβunsat/Acylsulfonamide | mTOR | +2.0 | Rapamycin-class allosteric FKBP12 binding |
| Adenosine Purine | Purine present | adenosine receptor | +0.5 | Purine is the defining adenosine scaffold |
| Kinase αβunsat warhead | α,β-unsat. carbonyl present | kinase | +0.5 | Covalent Michael acceptor warhead (EGFR) |
| Kinase sulfonamide-amine | Sulfonamide + TertAmine | kinase | +2.0 | Kinase linker hijacked by CA (Sulfonamide mw=2.0) |
| MAO warhead | (Propargylamine OR Hydrazine), no Sulfonamide/Nitrile/αβunsat | MAO | +2.5 | Irreversible MAO inhibitor: propargylamine→FAD adduct (selegiline/clorgiline) or hydrazine (phenelzine). Exclusions prevent CA/covalent-kinase false positives. Markers are routing-only (`_WARHEAD_ANNOTATIONS`, not in fg_database → no IDF shift) |
| Cysteine protease nitrile | Nitrile + Amide, no αβunsat/Pyrimidine/Quinazoline/Pyrrolopyrimidine/Fused-azolo/Sulfonamide | cysteine protease | +2.5 | Peptidomimetic nitrile warhead → thioimidate with catalytic Cys25 (odanacatib/cathepsin-K class). Exclusions strip kinase-hinge / covalent-kinase / CA contexts. In 320-cpd benchmark matches exactly 12 compounds, all cysteine protease (zero collision) |

\### Pyrimidine router (`_pyrimidine_router`, utils/target_predictor.py)

Mutually-exclusive routing for diazine-bearing ATP-pocket / purine-mimetic binders.
Replaces the old standalone morpholino-diazine mTOR rule. Evaluated in order; **exactly one branch fires** (conflict-free by construction). Gated on Pyrimidine OR Triazine present.

| Branch | Condition (after no earlier branch fired) | Target | Bonus | Benchmark provenance |
|---|---|---|---|---|
| 1 | Morpholine present | mTOR | +2.0 | 14 compounds, all mTOR. Morpholine O H-bonds hinge Val2240 |
| 2 | Fused azolo-diazine core (purine-mimetic) | adenosine receptor | +2.0 | 13 compounds = 12 adenosine + 1 mTOR (already-miss). Also suppresses CYP450 (negative constraint) |
| 3 | mono-Pyrimidine, no Methylsulfone/Hydroxamate/COOH/Aldehyde/Steroid | kinase | +2.0 | kinase-dominant; exclusions protect COX(Methylsulfone)/HDAC(Hydroxamate)/GPCR(COOH/Aldehyde)/NR(Steroid) HITs |

Branch-3 exclusions are competing pharmacophores whose own FG votes/rules already claim the compound. Branch 3 also short-circuits when `_has_free_heme_azole` is true, so a free-triazole antifungal with a *separate* fluoropyrimidine (e.g. **voriconazole**) falls through to the CYP450 azole rule instead of being misrouted to kinase (verified: voriconazole & fluconazole → CYP450; zero benchmark regression).

**Negative constraints** (suppress cytochrome P450 entirely):
- Hydroxamate or Thiol present → Zn-chelation → HDAC/metalloprotease context
- Acylsulfonamide present → tubulin macrolide warhead context
- **Fused azolo-diazine + Pyrimidine** → purine-mimetic core; ring N locked in fused diazine cannot coordinate heme Fe (no CYP450 TP has this core)

---

\## Known structural limitations (do NOT try to fix with mw tuning)

1. **mTOR 85% (17/20)**: SIROLIMUS fixed by macrolide rule; 16 ATP-competitive TORKinibs fixed by morpholino-diazine rule (2026-06-15). 3 remaining have NO morpholine: SAPANISERTIB & CHEMBL3645910 (pyrimidine core only, no morpholine), CHEMBL3681183 (Hydroxyl+Imidazole → CYP450). These need fused-N-heteroaromatic core detection (方案 4, flagged in Next tasks) — do NOT try to fix with the morpholino rule.
2. **Adenosine receptor 60% (12/20)**: +7 by pyrimidine router branch 2 (fused-azolo-diazine→adenosine, 2026-06-15). 8 remaining have NO purine-mimetic fused core: sparse Phenol+Phenyl+Halogen→NR (CHEMBL2024114, 97760), Nitrile+Triazole→cys protease/kinase (CHEMBL5177144, 5171044), Steroid (CHEMBL369573), Thiazole+Nitrile→CYP (CHEMBL2419137, 2419150), non-fused Thiazole+Pyrimidine (CHEMBL3917647).
3. **CYP450 95% (19/20)**: 1 remaining failure = TAZAROTENIC ACID (pyridine scaffold; no azole/halogen FG → invisible to CYP scoring). Thiazole SMARTS (2026-06-04) fixed ritonavir-class ×5 compound.
4. **Serine protease 65% (13/20)**: +1 by Guanidine FG (CHEMBL353760, 2026-06-16). 7 remaining peptidomimetics have NO S1 Arg-mimetic at all (no Benzamidine/Guanidine/amidine) — verified; structural, not a missing-pattern gap. They look like NR/tubulin/GPCR/CA.
5. **Kinase 90% (18/20)**: +4 by pyrimidine router branch 3 (mono-pyrimidine→kinase, 2026-06-15; recovered ERLOTINIB, CHEMBL29197/176582/174426). 2 remaining: CHEMBL5270693 (strong GPCR score 6.31 from TertAmine+Indole+Phenyl beats kinase bonus), CHEMBL4537790 (Steroid scaffold → androgen).
6. **NR 80% (16/20)**: 2 Acylsulfonamide→tubulin (irreconcilable without hurting tubulin). 2 purely structural.

---

\## Confirmed negative results (do NOT retry)

- **Halogen mw boost (1.2/1.5/2.0)**: All tested 2026-06-01. mw=1.2 → net 0; mw=1.5/2.0 → net -1 (gains 2 CYP450, loses 3: kinase+NR+HDAC). **Halogen is a promiscuity feature, not CYP450-specific.**
- **Macrolide mw=1.2**: net 0. SIROLIMUS Ketone+Lactone HDAC score always beats Macrolide mTOR score.
- **GPCR saturation (diminishing returns)**: Risky — GPCR compounds with Indole+Phenol would lose to serotonin receptor (IDF=3.555) when GPCR second vote is reduced below 3.555. Current 100% GPCR accuracy depends on full FG accumulation.
- **Carboxylic acid → NR conditional**: Adding CA+acid → NR would break INDOMETHACIN (COX HIT) which has identical FG profile (CA+acid+Ether+Phenyl+Halogen) to the NR compound CHEMBL2323507.
- **Stable sort for NR/tubulin tie-breaking**: Using `kind='stable'` in sort_values swaps GS-9256 (tubulin HIT → MISS) for CHEMBL180681 (NR MISS → HIT). Net Top-1 = 0. Tubulin drops, which is worse.
- **Halogen → NR annotation**: Would fix CHEMBL180681 (Phenyl+Halogen→NR tie) but breaks INDOMETHACIN COX HIT (Halogen+Phenyl+Ether+COOH → NR would beat COX). Net negative.
- **Benzimidazole with kinase/tubulin target classes**: Adding "kinase"+"tubulin" to Benzimidazole target classes (2026-06-04 attempt) decreased kinase IDF (5→6 annotators) and tubulin IDF (4→5 annotators), flipping 6 kinase→GPCR and 7 tubulin→CA. Net -13. Fix: Benzimidazole target_classes=[] (scaffold marker only). **Do NOT add kinase or tubulin to Benzimidazole annotations.**
- **GS-9256 (tubulin compound with thiazole ring)**: Cannot be distinguished from ritonavir-class at FG level (both have Thiazole+Phenyl+Ether profile). Thiazole addition causes GS-9256 to be incorrectly predicted as CYP450 (net -1 tubulin). Accept as structural limitation.

---

\## Next tasks (prioritised by ROI)

\### Blind-spot classes (7 classes with no FG rules) — phased plan

**Tier 1 DONE (2026-06-15):** MAO + COMT added to benchmark (260 cpds, 13 classes).
- COMT 8/20 (nitrocatechol via existing Phenol+Catechol — already maxed; rest are research series).
- MAO 2/20 (propargylamine/hydrazine warhead rule). pChEMBL bias caps both; see per-class note.

**Tier 2 PARTIAL (2026-06-15):** benchmark extended to 16 classes / 320 cpds (compounds.csv,
gitignored). Identities verified via unified map.
- **cysteine protease DONE: 12/20** — Nitrile+Amide gated rule (`_cysteine_protease_conditional_bonus`),
  zero collision, zero regression on core 11 (188/220). odanacatib/cathepsin-K nitrile-warhead class.
- **topoisomerase DONE: 5/20** — Anthraquinone added as a real voting FG (FG_SMARTS + fg_database,
  mw=2.5, known_target_classes=[topoisomerase]; table rebuilt). As topo's sole annotator its IDF≈3.7,
  so anthracyclines (doxorubicin etc.) beat the Ketone→HDAC pull. Anthraquinone matches exactly 5
  benchmark compounds, all topo (zero collision). N 41→42 IDF shift caused zero regression.
- **xanthine oxidase 0/20 — skip**: pChEMBL-bias, research series (Phenol+Amide+Pyrimidine), no
  allopurinol/febuxostat; pyrimidine→kinase collision. Structural limit.

**Tier 3 (deprioritised):** PDE (too generic + collisions), ribosome (only 3 compounds + hard
aminoglycosides). Likely structural limits like adenosine/serine protease.

\### 3D-fallback layer (skeleton DONE, real models pending)

`utils/target_predictor.py`: `assess_confidence()` + `register_fallback_3d()` +
`_finalize_with_fallback()`. Predictions get a `confidence` tag (high/low/none) from
observable signal only (top1 score < 3.0 OR top1−top2 margin < 0.75 → low; empty votes → none).
Low/none route to a pluggable 3D hook; default `_stub_fallback_3d` returns None → **zero
regression** (verified: core 11 stays 188/220). Plug real models via `register_fallback_3d`.

**Gate-selectivity finding (320-cpd benchmark)** — the confidence gate is SAFE but PARTIAL:
- Routes 42/105 misses (40%); also flags 41/215 hits (would be at risk only if a real model overrides).
- **63/105 misses are HIGH confidence (FG confidently wrong) → NOT reachable by the gate**:
  NO-ANSWER 39 (e.g. XO→kinase via pyrimidine), FALSE-POS 24 (e.g. COMT loses to CA/NR at score 6+).
- Implication: a confidence gate alone cannot trigger the FALSE-POSITIVE re-rank cases (the wrong
  winner is *strong*). Reaching the 63 confident-wrong misses needs a parallel 3D path + meta-reconciler
  that may override high-confidence FG calls — which carries regression risk and must be validated hard.
  Method map (from failure-mode analysis): shape/ProLIF = retrieval for NO-ANSWER; Gnina rescore =
  re-rank for FALSE-POS (but needs an always-on or top-K trigger, not the confidence gate).

**ProLIF fallback IMPLEMENTED (`utils/fallback_3d.py`)**: `Fallback3D` base + `ProLIFFallback` with the
full 4-step pipeline written: `_embed_3d` (RDKit ETKDGv3+MMFF) → `_dock` (smina subprocess,
`--autobox_ligand` on the reference co-crystal ligand) → `_compute_ifp` (ProLIF Fingerprint vs receptor)
→ `_match_reference` (Jaccard on interaction-key sets). `build_override()` re-rank/merge contract tested.
**IFP representation = set of "PROTRES.Interaction" keys** (e.g. ASP189.HBDonor) shared by builder and
fallback — NOT raw bitvectors (those aren't aligned across molecules); Tanimoto = Jaccard on key-sets.
`propose()` swallows all errors → None, and short-circuits when no reference library exists, so it stays
**zero-regression** until a library is built (verified 0/320 top1 change registered). Heavy deps
(ProLIF/MDAnalysis/docking) lazy-imported. **Runs once installed + reference built.** Regression surface
for Option 1 = the 41 low-confidence HITS (overrides must not break them).

**Reference-library builder DONE (`utils/build_prolif_reference.py check|build`)**: builds
db/prolif_reference_ifp.json from PDB **co-crystals** (real poses → ProLIF IFP, no docking on the
reference side; docking only at query time). Seed set = serine-protease complexes (trypsin 3PTB+BEN,
thrombin 1OYT/1DWD, factor Xa 2ZFF/1F0R — expand with more non-benzamidine peptidomimetics). Heavy
deps lazy; `check` runs anywhere. **Env status (2026-06-16): chem_target has rdkit + ProLIF 2.1.0 +
MDAnalysis 2.9.0 + requests ✓ (build-ready); no docking backend yet (query-time, `conda install -c
conda-forge smina`).** NOTE: the user's interactive shell can't `conda activate` (run via the full
path `C:\Users\User\miniconda3\envs\chem_target\python.exe`). db/prolif_pdb/ + db/prolif_reference_ifp.json
gitignored.

**⛔ Build blocker (protein prep) — reference library NOT yet built (2026-06-16):** the builder now
auto-fixes conda OpenBabel's plugin path (`_setup_openbabel`: BABEL_LIBDIR/DATADIR + add_dll_directory
from sys.prefix), protonates ligand + protein, and sources them separately (OpenBabel drops HETATM
resnames). BUT OpenBabel mangles **protein** chain/residue topology → ProLIF raises
`KeyError: ResidueId(...)` in `Fingerprint.generate`. **Fix: prepare the protein with a
topology-preserving tool (PDBFixer/openmm or `reduce`), not OpenBabel** — `conda install -c conda-forge
pdbfixer openmm`, then swap the protein branch of `_ifp_from_complex` to a PDBFixer addMissingHydrogens
pass (keep OpenBabel only for the ligand). Decision pending: this is real structural-prep yak-shaving
for a ~7-compound (serine-protease) payoff; the FG core (188/220) is unaffected and all 3D-fallback
infra is committed, so the ProLIF PoC can resume anytime without blocking other work.

**3D-fallback tune list (once installed + library built):**
- **Receptor prep**: `_dock` currently feeds smina a raw PDB receptor (no PDBQT prep). Fine for a PoC,
  but for docking accuracy add proper prep (meeko / `prepare_receptor` → PDBQT: protonation, charges,
  rigid/flex split). Same for the autobox ligand reference.
- **Confidence gate thresholds** (`target_predictor`): `_CONF_MIN_TOP1_SCORE=3.0`, `_CONF_MIN_MARGIN=0.75`
  — currently routes 40% of misses + flags 41 hits; tune vs the hit/miss cross-tab.
- **ProLIF similarity→score** (`ProLIFFallback`): `sim_threshold=0.6`, `score_scale=8.0` — sim must be
  high to be competitive (FALSE-POS recovery needs sim≈0.9+ to beat score-6 winners; NO-ANSWER ≈0.7).
- **Reference set coverage**: seed has 5 serine-protease co-crystals; add more non-benzamidine
  peptidomimetic complexes so the 8 SP misses have a near neighbour.
- **Docking determinism**: smina `--seed 0` already set; consider `--num_modes>1` + best-IFP-over-poses,
  and tune `--exhaustiveness`/`--autobox_add`.
- **Override policy**: `build_override` currently boosts/inserts one proposal; for FALSE-POS re-rank at
  scale, consider top-K docking + a meta-reconciler (Option 2) rather than the single low-conf gate.

\### Other improvements

1. ~~SDF / MOL2 / SMI / InChI input~~ ✅ DONE (2026-06-16, io_handler.py: RDKit for SDF/SMI/InChI, OpenBabel for MOL2)
2. ~~Shape descriptors~~ ✅ DONE (2026-06-16) as an **informational** layer (`utils/shape_descriptors.py`,
   opt-in `predict --shape` → NPR1/NPR2, Rg, asphericity, eccentricity, spherocity, shape_class).
   **⚠ The original hypothesis was BACKWARDS** — benchmark medians show **GPCR** ligands are the LARGE,
   elongated ones (Rg 8.5, asphericity 0.73), while **CYP450** is compact/globular (Rg 4.7, asph 0.36).
   Shape does carry class signal (GPCR largest; NR/COMT smallest Rg≈3.6; XO/adenosine flat asph≈0.7),
   but the most shape-separable class (GPCR) is already 100%, so **no shape SCORING rule was added**
   (would risk the 188 core for no gain). Shape-based help for the misses (e.g. flat XO vs kinase)
   would need the same gated 3D-fallback treatment as ProLIF — deferred.
3. ~~Serine protease guanidine coverage~~ ✅ DONE (2026-06-16): Guanidine FG (mw 3.0→SP) recovered
   CHEMBL353760 (SP 12→13); zero regression (the 4 guanidine compounds were all already misses), and
   the IDF shift incidentally nudged NR 16→17. Data finding: only 1/8 SP misses had guanidine; the
   other 7 carry NO S1 Arg-mimetic → structural limit, no further pattern will help.
4. **ProLIF 3D-fallback (PAUSED)** — infra committed; resume needs PDBFixer protein-prep (see 3D-fallback section)

\### ✅ 方案 4 (DONE — detection complete; no further scoring ROI on this benchmark)

**Functional core (2026-06-15)**: `_detect_fused_azolo_diazine` (aromatic 5-ring≥2N fused to
6-ring≥2N = purine / triazolopyrimidine / pyrazolopyrimidine). Drives the pyrimidine router branch 2
(→adenosine, +7) and the CYP450 negative constraint. Recipe that worked: **routing-only Python
detector, NOT in fg_database.json → zero IDF impact** (sidesteps the Benzimidazole −13 lesson).

**Annotation cores (2026-06-15)**: `_SCAFFOLD_ANNOTATIONS` in fg_detector.py adds Quinazoline /
Pyrrolopyrimidine / Pyridopyrimidine / Benzoxazole as **annotation-only** labels (not in FG_SMARTS,
not in fg_database, not consumed by any rule). Purpose: scaffold cores no longer show up merely as
"Phenyl ring" in reports; infrastructure for future rules. **Zero benchmark change** (cast no votes).

**ROI finding (verified 2026-06-15)** — do NOT add scoring rules for remaining cores:
Of the 32 misses at 188/220, only 4 carry an undetected fused-N core and **none are fixable by core
detection**: CHEMBL4108739/353760 (serine protease peptidomimetics — benzoxazole/benzothiazole, but
no SP signal), SAPANISERTIB (mTOR — benzoxazole; routing benzoxazole→mTOR is promiscuous, fixes only
1), GS-9256 (tubulin — quinoline; documented irreconcilable). The remaining misses are structural
(no Benzamidine, Acylsulfonamide→tubulin, sparse Phenol+Halogen, Steroid), not core-detection-limited.
The azolo-diazine subset already captured all fused-core scoring ROI this benchmark offers.

**Also**: 方案 3 standalone-pyrimidine-voting is SUPERSEDED by the router (branch 3 does it safely
with exclusions). Do not add Pyrimidine to fg_database known_target_classes (would disrupt IDF).

