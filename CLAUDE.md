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
| v1/v2 prototype | ✅ Done |
| Layered refactor (constants/db/utils/main) | ✅ Done |
| 35 FG SMARTS (`constants/fg_smarts.py`) | ✅ Done |
| `db/fg_database.json` (35 entries, mechanistic_weight) | ✅ Done |
| `db/fg_residue_table.csv` (BioLiP rebuild) | ✅ Done |
| `db/residue_3d_poses.json` + `db/local_env/*.sdf` | ✅ Done |
| `utils/target_predictor.py` (IDF × mechanistic_weight) | ✅ Done |
| `utils/report_generator.py` (HTML individual + batch) | ✅ Done |
| `run_benchmark.py` (11-class × 20-compound curated) | ✅ Done |
| **Benchmark Top-1: 137/220 = 62.3%** | ✅ Current best |
| **Benchmark Top-3: 148/220 = 67.3%** | ✅ Current best |
| SDF / MOL2 input support | 🔲 Pending |
| CYP450 conditional motif scoring | 🔲 Next |
| Negative constraint rules (kinase/HDAC suppression) | 🔲 Next |
| Shape / physicochemical descriptors | 🔲 Future |
| Merge dev/validation → master | 🔲 Pending |

---

\## Benchmark per-class results (current best, branch: dev/validation)

| Class | Top-1 | Top-3 | Notes |
|---|---|---|---|
| GPCR | 20/20 = 100% | 20/20 | ✅ |
| HDAC | 20/20 = 100% | 20/20 | ✅ Fixed by Ketone mw=2.0 |
| Carbonic anhydrase | 20/20 = 100% | 20/20 | ✅ |
| Tubulin | 20/20 = 100% | 20/20 | ✅ Fixed by Acylsulfonamide mw=2.0 |
| Nuclear receptor | 17/20 = 85% | 20/20 | 3 losses: Acylsulfonamide→tubulin conflict |
| Serine protease | 12/20 = 60% | 12/20 | 8 failures: no Benzamidine FG signal |
| COX | 11/20 = 55% | 12/20 | 9 failures: Sulfonamide→CA hijack |
| Kinase | 8/20 = 40% | 12/20 | Structural limit: ATP-pocket overlap |
| CYP450 | 7/20 = 35% | 8/20 | Structural limit: low IDF, weak FG signal |
| Adenosine receptor | 4/20 = 20% | 5/20 | Structural limit: 16/20 non-purine |
| mTOR | 0/20 = 0% | 0/20 | Structural limit: 19/20 ATP-competitive |

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

\## Known structural limitations (do NOT try to fix with mw tuning)

1. **mTOR 0%**: SIROLIMUS Ketone+Lactone → HDAC score (5.838) > Macrolide→mTOR (4.266). Irreconcilable without conditional scoring.
2. **Adenosine receptor 20%**: 16/20 compounds lack Purine/Imidazole/Xanthine scaffold.
3. **Kinase 40%**: ATP-pocket FGs (Phenyl/Amide/Ether) shared with GPCR/NR/CYP450.
4. **CYP450 35%**: IDF too low (n=7, IDF=1.609). Global Halogen mw boost tested and **rejected** (+2 CYP450 but -3 other classes, ROI negative).

---

\## Confirmed negative results (do NOT retry)

- **Halogen mw boost (1.2/1.5/2.0)**: All tested 2026-06-01. mw=1.2 → net 0; mw=1.5/2.0 → net -1 (gains 2 CYP450, loses 3: kinase+NR+HDAC). **Halogen is a promiscuity feature, not CYP450-specific.**
- **Macrolide mw=1.2**: net 0. SIROLIMUS Ketone+Lactone HDAC score always beats Macrolide mTOR score.

---

\## Next tasks (prioritised by ROI, from onboard/cyp450.md)

\### Highest ROI — implement first

1. **CYP450 conditional motif scoring** (`onboard/cyp450.md` Phase 2)
   - Replace single-FG Halogen/Imidazole scoring with combo rules:
     - `Imidazole + hydrophobic tail` → CYP450
     - `Halogenated phenyl + oxidation site` → CYP450
   - Implementation: new `_cyp450_conditional_score()` in `utils/target_predictor.py`

2. **Negative constraint rules** (`onboard/cyp450.md` Phase 5)
   - `flat hinge geometry + dual HBD/HBA + rigid aromatic` → suppress CYP450 (kinase filter)
   - `Zn chelation (Hydroxamate/Thiol)` → suppress CYP450 (HDAC filter)

3. **GPCR saturation / diminishing return** (`onboard/cyp450.md` Phase 6 Step 11)
   - First amine/aromatic match = 1.0, second = 0.3, third = 0.1
   - Prevents GPCR from overwhelming mTOR/kinase with repeated generic FGs

\### Medium ROI

4. **Shape descriptors** (Phase 3): PMI, radius of gyration → CYP450 = elongated ligand
5. **Oxidation soft spot detector** (Phase 4): benzylic C, allylic position, terminal methyl

\### Admin / housekeeping

6. Merge `dev/validation` → `master` once CYP450 conditional scoring is stable
7. Update `db/fg_residue_table.csv` if new FGs are added
8. SDF / MOL2 input support (`utils/io_handler.py`)

---

\## Immediate tasks (in order)

1\. Implement CYP450 conditional motif in `utils/target_predictor.py`
2\. Implement negative constraints (kinase + HDAC suppression of CYP450)
3\. Re-run `run_benchmark.py run` and compare against 137/220 baseline
4\. If net positive: commit + push on `dev/validation`
5\. Merge to `master` when satisfied

