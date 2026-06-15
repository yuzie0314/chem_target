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
│   └── fg_smarts.py    # FG_SMARTS: name → SMARTS (40 patterns + Steroid Python = 41 total)
├── db/             # Auto-generated. Never hand-edit.
│   ├── fg_database.json       # 41 FG metadata: smarts/targets/mechanistic_weight
│   ├── fg_residue_table.csv   # 40 SMARTS + Steroid × 20 AA BioLiP co-occurrence matrix
│   ├── ccd_smiles_cache.json  # RCSB CCD SMILES cache
│   └── residue_3d_poses.json  # Cα + ligand centroid 3D records
├── utils/          # Pure functions. No side effects where possible.
│   ├── fg_detector.py         # detect_smarts(), _detect_steroid_core(), _detect_fused_azolo_diazine()
│   ├── target_predictor.py    # IDF × mw scoring + conditional rules + _pyrimidine_router()
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



Currently supported:

\- CSV: col1 = compound name, col2 = SMILES



Planned (extend io\_handler.py):

\- SDF files

\- MOL2 files

\- InChI strings



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
| 40 FG SMARTS + Steroid Python = 41 total (`constants/fg_smarts.py`) | ✅ Done |
| `db/fg_database.json` (41 entries incl. Triazole+Thiazole+Benzimidazole+Morpholine+Pyrimidine+Triazine, mechanistic_weight) | ✅ Done |
| `db/fg_residue_table.csv` (BioLiP rebuild, 40 SMARTS + Steroid columns) | ✅ Done |
| `db/residue_3d_poses.json` + `db/local_env/*.sdf` | ✅ Done |
| `utils/target_predictor.py` (IDF × mechanistic_weight) | ✅ Done |
| `utils/report_generator.py` (HTML individual + batch) | ✅ Done |
| `run_benchmark.py` (11-class × 20-compound curated) | ✅ Done |
| **Core 11-class Top-1: 188/220 = 85.5%** (mechanistic classes) | ✅ Current best |
| **Core 11-class Top-3: 196/220 = 89.1%** | ✅ Current best |
| **Extended 13-class Top-1: 198/260 = 76.2%** (incl. MAO+COMT) | ✅ |
| MAO covalent-warhead rule (Propargylamine/Hydrazine) | ✅ Done |
| COMT (nitrocatechol via existing Phenol+Catechol) | ✅ 8/20 (pChEMBL-bias limited) |
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
| 方案 4 fused-N core: azolo-diazine detector (functional, +7) + Quinazoline/Pyrrolopyrimidine/Pyridopyrimidine/Benzoxazole (annotation-only) | ✅ Done |
| SDF / MOL2 input support | 🔲 Pending |
| Shape / physicochemical descriptors | 🔲 Future |

---

\## Benchmark per-class results (current best, branch: dev/validation)

| Class | Top-1 | Top-3 | Notes |
|---|---|---|---|
| GPCR | 20/20 = 100% | 20/20 | ✅ |
| HDAC | 20/20 = 100% | 20/20 | ✅ |
| Carbonic anhydrase | 20/20 = 100% | 20/20 | ✅ |
| Tubulin | 19/20 = 95% | 19/20 | -1 GS-9256 (thiazole+ether → CYP rule; profile indistinguishable from ritonavir-class) |
| Nuclear receptor | 16/20 = 80% | 20/20 | 4 losses: 2× Acylsulfonamide→tubulin + 2× structural |
| Serine protease | 12/20 = 60% | 12/20 | 8 failures: no Benzamidine FG signal |
| COX | 15/20 = 75% | 17/20 | Fixed +4 by Indole+Sulfonamide motif |
| Kinase | 18/20 = 90% | 18/20 | +4 by pyrimidine router (mono-pyrimidine→kinase, branch 3); earlier +6 αβunsat+Sulfonamide. 2 remaining: 1 strong-GPCR (CHEMBL5270693), 1 Steroid |
| CYP450 | 19/20 = 95% | 19/20 | Fixed +12 total; 5 ritonavir-class by Thiazole SMARTS; 1 TAZAROTENIC ACID structural. Pyrimidine guard added (no CYP450 TP has pyrimidine) |
| Adenosine receptor | 12/20 = 60% | 12/20 | +7 by pyrimidine router (fused-azolo-diazine→adenosine, branch 2). 8 remaining: no purine-mimetic core (Phenol/Halogen sparse, or Nitrile/Steroid) |
| mTOR | 17/20 = 85% | 17/20 | Fixed +16 by morpholino-diazine motif (16 ATP-competitive TORKinibs); +SIROLIMUS by macrolide rule. 3 remaining have no morpholine (SAPANISERTIB, CHEMBL3645910, CHEMBL3681183) |
| COMT | 8/20 = 40% | — | nitrocatechol (entacapone/opicapone) via existing Phenol+Catechol; other 12 = research series w/o nitrocatechol (pChEMBL-bias) |
| MAO | 2/20 = 10% | — | propargylamine (clorgiline) via MAO warhead rule; 18 = single research series (Sec/Tert amine, no MAO pharmacophore) |

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
| Hydroxamate | 2.5 | HDAC | Bidentate Zn chelation |
| Sulfonamide | 2.0 | carbonic anhydrase | Zn coordination |
| Acylsulfonamide | 2.0 | tubulin | Macrolide warhead |
| Ketone | 2.0 | HDAC | α-keto warhead in HDAC natural products |
| Steroid | 2.0 | nuclear receptor (+ subtypes) | Steroidal scaffold → NR |
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
4. **Serine protease 60% (12/20)**: 8 failures have no Benzamidine. Peptidomimetics look like NR/tubulin/GPCR.
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

**Tier 2 (next, if pursued):** topoisomerase (Anthraquinone + quinolone-3-COOH SMARTS),
xanthine oxidase (Pyrazolopyrimidine already detected + Thiazole-COOH), cysteine protease
(Nitrile warhead — but collides with kinase, must gate). The ChEMBL target maps are now unified
(done 2026-06-15); still verify each class's downloaded compounds are the right drug class (pChEMBL
bias often returns research analogs without the textbook pharmacophore) before adding rules.

**Tier 3 (deprioritised):** PDE (too generic + collisions), ribosome (only 3 compounds + hard
aminoglycosides). Likely structural limits like adenosine/serine protease.

\### Other improvements

1. **SDF / MOL2 input support** (`utils/io_handler.py`) — currently only CSV supported
2. **Shape descriptors** (PMI, radius of gyration) — would help distinguish CYP450 elongated ligands from compact GPCR ligands
3. **Serine protease Benzamidine coverage** — 8 failures have no Benzamidine (peptidomimetics); possible solution: add guanidine or charged amidino group pattern

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

