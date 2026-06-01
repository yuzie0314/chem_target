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

chem\_target/

├── constants/      # Static lookup tables only. No logic.

│   ├── fg\_names.py     # FG\_NAMES: rdkit code → human readable name

│   └── fg\_smarts.py    # FG\_SMARTS: name → SMARTS pattern

├── db/             # Auto-generated. Never hand-edit.

│   └── fg\_database.json    # Fetched from PubChem + ChEMBL

├── utils/          # Pure functions. No side effects where possible.

│   ├── io\_handler.py       # File reading (CSV, SDF, etc.)

│   ├── fg\_detector.py      # Functional group detection logic

│   ├── db\_updater.py       # Fetch + update fg\_database.json

│   └── visualizer.py       # RDKit drawing + SVG output

├── data/           # User input files

├── output/         # Generated output — CSV tables, images, reports

└── main.py         # Entry point only. Thin. No logic here.

```



\*\*Never put constants inside logic files.\*\*

\*\*Never put logic inside main.py.\*\*

\*\*Never hardcode values that belong in constants/ or db/.\*\*



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

\- Do not write to db/ manually — always via db\_updater.py

\- Do not over-engineer early — get it working first, then clean up

\- Do not ignore OpenBabel for format handling — it is already installed



\---



\## Current status (2026-06-01)

| Component | Status |
|---|---|
| conda env + rdkit + openbabel | ✅ Done |
| 35 FG SMARTS (`constants/fg_smarts.py`) | ✅ Done |
| `db/fg_database.json` (35 entries, mechanistic_weight) | ✅ Done |
| `db/fg_residue_table.csv` (BioLiP rebuild) | ✅ Done |
| `db/residue_3d_poses.json` + `db/local_env/*.sdf` | ✅ Done |
| `utils/target_predictor.py` (IDF × mechanistic_weight) | ✅ Done |
| `utils/report_generator.py` (HTML individual + batch) | ✅ Done |
| `run_benchmark.py` (11-class × 20-compound curated) | ✅ Done |
| **Benchmark Top-1: 149/220 = 67.7%** | ✅ Current best |
| **Benchmark Top-3: 157/220 = 71.4%** | ✅ Current best |
| CYP450 conditional motif scoring (azole rule) | ✅ Done |
| Negative constraint rules (Hydroxamate/Thiol/Acylsulfonamide → suppress CYP450) | ✅ Done |
| COX indole-sulfonamide motif | ✅ Done |
| mTOR macrolide conditional motif | ✅ Done |
| Adenosine receptor Purine bonus | ✅ Done |
| Kinase α,β-unsat carbonyl covalent warhead bonus | ✅ Done |
| SDF / MOL2 input support | 🔲 Pending |
| Shape / physicochemical descriptors | 🔲 Future |
| Merge dev/validation → master | 🔲 Pending |

---

\## Benchmark per-class results (current best, branch: dev/validation)

| Class | Top-1 | Top-3 | Notes |
|---|---|---|---|
| GPCR | 20/20 = 100% | 20/20 | ✅ |
| HDAC | 20/20 = 100% | 20/20 | ✅ |
| Carbonic anhydrase | 20/20 = 100% | 20/20 | ✅ |
| Tubulin | 20/20 = 100% | 20/20 | ✅ |
| Nuclear receptor | 16/20 = 80% | 20/20 | 4 losses: 2× Acylsulfonamide→tubulin + 2× structural |
| Serine protease | 12/20 = 60% | 12/20 | 8 failures: no Benzamidine FG signal |
| COX | 15/20 = 75% | 17/20 | Fixed +4 by Indole+Sulfonamide motif |
| Kinase | 13/20 = 65% | 14/20 | Fixed +5 by α,β-unsat covalent warhead motif |
| CYP450 | 7/20 = 35% | 8/20 | Structural limit: 13 failures |
| Adenosine receptor | 5/20 = 25% | 5/20 | Fixed +1 by Purine bonus; 15 structural |
| mTOR | 1/20 = 5% | 1/20 | Fixed SIROLIMUS; 19 ATP-competitive structural |

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
| All others | 1.0 | — | Default |

---

\## Conditional scoring rules (utils/target_predictor.py)

All rules are pre-IDF bonuses (multiplied by IDF before adding to final score).

| Rule | Condition | Target | Bonus | Rationale |
|---|---|---|---|---|
| CYP450 azole | Imidazole + {Phenyl/Ether/Halogen}, no Ketone/Purine/αβunsat | cytochrome P450 | +2.0 | Azole antifungal heme-Fe coordination |
| COX indole-sulfonamide | Indole + Sulfonamide | COX | +2.0 | Indole scaffold + COX-2 selectivity pocket |
| mTOR macrolide | Macrolide, no Thiol/αβunsat/Acylsulfonamide | mTOR | +2.0 | Rapamycin-class allosteric FKBP12 binding |
| Adenosine Purine | Purine present | adenosine receptor | +0.5 | Purine is the defining adenosine scaffold |
| Kinase warhead | α,β-unsat. carbonyl present | kinase | +0.5 | Covalent Michael acceptor warhead (EGFR etc.) |

**Negative constraints** (suppress cytochrome P450 entirely):
- Hydroxamate or Thiol present → Zn-chelation → HDAC/metalloprotease context
- Acylsulfonamide present → tubulin macrolide warhead context

---

\## Known structural limitations (do NOT try to fix with mw tuning)

1. **mTOR 5% (1/20)**: Only SIROLIMUS (macrolide) fixed. 19/20 ATP-competitive inhibitors look like kinase/NR compounds with no mTOR-specific FG signal.
2. **Adenosine receptor 25% (5/20)**: 15/20 failures lack Purine/Xanthine scaffold entirely; structurally indistinguishable from NR/tubulin/cysteine protease compounds.
3. **CYP450 35% (7/20)**: 13 failures. 4 CYP450 compounds are aryl-COOH drugs (same FG profile as NSAIDs, predict COX). Triazole-class azoles (fluconazole, voriconazole) cannot be fixed without adding Triazole SMARTS.
4. **Serine protease 60% (12/20)**: 8 failures have no Benzamidine. These peptidomimetics (factor Xa, thrombin inhibitors) look like GPCR/NR/tubulin compounds.
5. **Kinase 65% (13/20)**: 7 remaining failures. 2 stolen by CYP450 azole rule (Imidazole+Phenyl+Halogen, no distinguishing FG). 1 has steroidal scaffold (androgen wins). 4 are sparse.

---

\## Confirmed negative results (do NOT retry)

- **Halogen mw boost (1.2/1.5/2.0)**: All tested 2026-06-01. mw=1.2 → net 0; mw=1.5/2.0 → net -1 (gains 2 CYP450, loses 3: kinase+NR+HDAC). **Halogen is a promiscuity feature, not CYP450-specific.**
- **Macrolide mw=1.2**: net 0. SIROLIMUS Ketone+Lactone HDAC score always beats Macrolide mTOR score.
- **GPCR saturation (diminishing returns)**: Risky — GPCR compounds with Indole+Phenol would lose to serotonin receptor (IDF=3.555) when GPCR second vote is reduced below 3.555. Current 100% GPCR accuracy depends on full FG accumulation.
- **Carboxylic acid → NR conditional**: Adding CA+acid → NR would break INDOMETHACIN (COX HIT) which has identical FG profile (CA+acid+Ether+Phenyl+Halogen) to the NR compound CHEMBL2323507.

---

\## Next tasks (prioritised by ROI)

\### Potentially achievable

1. **Add Triazole SMARTS** (`constants/fg_smarts.py`): fluconazole/voriconazole class (triazole antifungals) can't trigger CYP450 azole rule currently. Adding `c1cn[nH,n]1` as "Triazole" FG would cover these. Requires rebuilding BioLiP table.
2. **Kinase Sulfonamide hijack** (CHEMBL5594833): TerAmine+Sulfonamide+Phenyl+Halogen → CA(5.724) steals this kinase compound. A conditional: Sulfonamide + TerAmine + no Indole → kinase bonus? Low ROI (+1 kinase).

\### Admin / housekeeping

3. Merge `dev/validation` → `master`
4. Update `db/fg_residue_table.csv` if new FGs are added
5. SDF / MOL2 input support (`utils/io_handler.py`)

