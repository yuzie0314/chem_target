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
│   ├── fg_detector.py         # detect_smarts(), _detect_steroid_core()
│   ├── target_predictor.py    # IDF × mw scoring + conditional rules
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
| **Benchmark Top-1: 177/220 = 80.5%** | ✅ Current best |
| **Benchmark Top-3: 185/220 = 84.1%** | ✅ Current best |
| CYP450 conditional motif scoring (azole rule, Thiazole added) | ✅ Done |
| Negative constraint rules (Hydroxamate/Thiol/Acylsulfonamide → suppress CYP450) | ✅ Done |
| COX indole-sulfonamide motif | ✅ Done |
| mTOR macrolide conditional motif (rapalog) | ✅ Done |
| mTOR morpholino-diazine motif (ATP-competitive TORKinib) | ✅ Done |
| Adenosine receptor Purine bonus | ✅ Done |
| Kinase α,β-unsat carbonyl covalent warhead bonus | ✅ Done |
| Thiazole SMARTS + BioLiP table rebuild | ✅ Done |
| Benzimidazole SMARTS (scaffold detection only, no target votes) | ✅ Done |
| Morpholine + Pyrimidine + Triazine SMARTS (scaffold markers, no target votes) | ✅ Done |
| SDF / MOL2 input support | 🔲 Pending |
| Fused-N-heteroaromatic core descriptor (方案 4 — see Next tasks) | 🔲 Flagged |
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
| Kinase | 14/20 = 70% | 15/20 | Fixed +6 by αβunsat warhead + Sulfonamide+TertAmine |
| CYP450 | 19/20 = 95% | 19/20 | Fixed +12 total; 5 ritonavir-class by Thiazole SMARTS; 1 TAZAROTENIC ACID structural |
| Adenosine receptor | 5/20 = 25% | 5/20 | Fixed +1 by Purine bonus; 15 structural |
| mTOR | 17/20 = 85% | 17/20 | Fixed +16 by morpholino-diazine motif (16 ATP-competitive TORKinibs); +SIROLIMUS by macrolide rule. 3 remaining have no morpholine (SAPANISERTIB, CHEMBL3645910, CHEMBL3681183) |

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
| Morpholine | 1.0 | (none — scaffold marker) | PI3K/mTOR hinge-binder; promiscuous alone (gefitinib). Voting via morpholino-diazine conditional rule |
| Pyrimidine | 1.0 | (none — scaffold marker) | 1,3-diazine hinge-anchor; promiscuous alone. Voting via morpholino-diazine conditional rule |
| Triazine | 1.0 | (none — scaffold marker) | 1,3,5-triazine hinge-anchor (gedatolisib class). Voting via morpholino-diazine conditional rule |
| All others | 1.0 | — | Default |

---

\## Conditional scoring rules (utils/target_predictor.py)

All rules are pre-IDF bonuses (multiplied by IDF before adding to final score).

| Rule | Condition | Target | Bonus | Rationale |
|---|---|---|---|---|
| CYP450 azole | (Imidazole OR Triazole OR Thiazole) + {Phenyl/Ether/Halogen}, Ketone only if no Amide/TertAmine, no Purine/αβunsat/Sulfonamide | cytochrome P450 | +2.0 | Azole/triazole/thiazole heme-Fe coordination (fluconazole/ritonavir class) |
| CYP450 aryl-COOH A | COOH + Phenyl + Halogen, no Amide, no Ether | cytochrome P450 | +1.5 | Minimal aryl-halide CYP substrate |
| CYP450 aryl-COOH B | COOH + Amide + Ether + Phenyl + Halogen | cytochrome P450 | +1.5 | Extended aryl-halide CYP substrate |
| CYP450 ether-amine | Ether + TertAmine + Phenyl + Halogen, no Lactone/Amide/Nitrile | cytochrome P450 | +1.5 | CYP3A4 scaffold (aprepitant-type) |
| CYP450 amide-halide | Amide + Phenyl + Halogen, no Sulfonamide/COOH/Imidazole/αβunsat/Ether | cytochrome P450 | +0.6 | Minimal amide-halide CYP substrate (raised from 0.5 to compensate IDF shift from Thiazole) |
| COX indole-sulfonamide | Indole + Sulfonamide | COX | +2.0 | Indole scaffold + COX-2 selectivity pocket |
| mTOR macrolide | Macrolide, no Thiol/αβunsat/Acylsulfonamide | mTOR | +2.0 | Rapamycin-class allosteric FKBP12 binding |
| mTOR morpholino-diazine | Morpholine + (Pyrimidine OR Triazine) | mTOR | +2.0 | ATP-competitive TORKinib: morpholine O H-bonds hinge Val2240. In curated benchmark this combo = 16 compounds, all mTOR (zero collision) |
| Adenosine Purine | Purine present | adenosine receptor | +0.5 | Purine is the defining adenosine scaffold |
| Kinase αβunsat warhead | α,β-unsat. carbonyl present | kinase | +0.5 | Covalent Michael acceptor warhead (EGFR) |
| Kinase sulfonamide-amine | Sulfonamide + TertAmine | kinase | +2.0 | Kinase linker hijacked by CA (Sulfonamide mw=2.0) |

**Negative constraints** (suppress cytochrome P450 entirely):
- Hydroxamate or Thiol present → Zn-chelation → HDAC/metalloprotease context
- Acylsulfonamide present → tubulin macrolide warhead context

---

\## Known structural limitations (do NOT try to fix with mw tuning)

1. **mTOR 85% (17/20)**: SIROLIMUS fixed by macrolide rule; 16 ATP-competitive TORKinibs fixed by morpholino-diazine rule (2026-06-15). 3 remaining have NO morpholine: SAPANISERTIB & CHEMBL3645910 (pyrimidine core only, no morpholine), CHEMBL3681183 (Hydroxyl+Imidazole → CYP450). These need fused-N-heteroaromatic core detection (方案 4, flagged in Next tasks) — do NOT try to fix with the morpholino rule.
2. **Adenosine receptor 25% (5/20)**: 15/20 failures have no Purine/Xanthine at all; generic Phenyl/Halogen → NR/tubulin.
3. **CYP450 95% (19/20)**: 1 remaining failure = TAZAROTENIC ACID (pyridine scaffold; no azole/halogen FG → invisible to CYP scoring). Thiazole SMARTS (2026-06-04) fixed ritonavir-class ×5 compound.
4. **Serine protease 60% (12/20)**: 8 failures have no Benzamidine. Peptidomimetics look like NR/tubulin/GPCR.
5. **Kinase 70% (14/20)**: 6 remaining. 2 stolen by CYP450 azole (Imidazole+Phenyl+Halogen, no distinguisher). 1 Steroid scaffold. 3 sparse.
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

\### Next potential improvements

1. **SDF / MOL2 input support** (`utils/io_handler.py`) — currently only CSV supported
2. **Shape descriptors** (PMI, radius of gyration) — would help distinguish CYP450 elongated ligands from compact GPCR ligands
3. **Serine protease Benzamidine coverage** — 8 failures have no Benzamidine (peptidomimetics); possible solution: add guanidine or charged amidino group pattern

\### 🚩 方案 4 (FLAGGED — long-term task): Fused-N-heteroaromatic core descriptor

**Problem**: the SMARTS detector is blind to fused N-heteroaromatic *cores*. E.g. SAPANISERTIB's
pyrazolo-pyrimidine / benzoxazole core is detected only as "Phenyl ring"; the whole kinase-hinge
scaffold is invisible. This is the root cause of the 3 remaining mTOR misses (compounds with a
heteroaromatic hinge-binder but no morpholine) and likely contributes to kinase / adenosine sparsity.

**Proposed solution**: add a generalised fused-N-bicyclic-aromatic scaffold detector (analogous to
the Python `_detect_steroid_core` approach, or a family of SMARTS for pyrazolopyrimidine,
pyrrolopyrimidine, pyridopyrimidine, etc.). Would help mTOR (+up to 3), and possibly kinase /
adenosine receptor.

**Effort/risk**: HIGH — broad scaffold patterns risk cross-class IDF disruption (cf. Benzimidazole
−13 lesson). Must be added as scaffold markers (target_classes=[]) + gated conditional rules, and
validated with the full 220-compound benchmark before keeping. Defer until after SDF/MOL2 input.

**Also flagged**: 方案 3 (standalone Pyrimidine/aminopyrimidine voting) was rejected for now —
pyrimidine alone is too promiscuous (kinase) and only covers +2; revisit only inside 方案 4.

