# CLAUDE.md — chem_target project

## Who you are working with
Yu-Jen Chang — pharmacist, computational chemist, data scientist.
Deep functional group intuition from pharmacy background.
Primary docking tools: iGEMDOCK, SiMMap (not AutoDock/Glide).
Current stack: Python 3.10, RDKit, OpenBabel, pandas, conda.

---

## Project goal
A compound-to-target prediction tool (reverse docking / target fishing).
Input: molecule (SMILES or SDF). Output: predicted enzyme/protein targets
based on functional group analysis. Long-term: consultancy tool + SaaS.

---

## Architecture rules — always follow these

```
chem_target/
├── constants/      # Static lookup tables only. No logic.
│   ├── fg_names.py     # FG_NAMES: rdkit code → human readable name
│   └── fg_smarts.py    # FG_SMARTS: name → SMARTS pattern
├── db/             # Auto-generated. Never hand-edit.
│   └── fg_database.json    # Fetched from PubChem + ChEMBL
├── utils/          # Pure functions. No side effects where possible.
│   ├── io_handler.py       # File reading (CSV, SDF, etc.)
│   ├── fg_detector.py      # Functional group detection logic
│   ├── db_updater.py       # Fetch + update fg_database.json
│   └── visualizer.py       # RDKit drawing + SVG output
├── data/           # User input files
├── output/         # Generated output — CSV tables, images, reports
└── main.py         # Entry point only. Thin. No logic here.
```

**Never put constants inside logic files.**
**Never put logic inside main.py.**
**Never hardcode values that belong in constants/ or db/.**

---

## Coding conventions

- Python 3.10
- Type hints on all function signatures
- Docstrings on all functions (one-line minimum)
- Constants in UPPER_SNAKE_CASE
- Functions in snake_case
- Early returns over nested if-else
- Explicit is better than implicit — no magic numbers

---

## DB update strategy

- Primary source: PubChem API (`https://pubchem.ncbi.nlm.nih.gov/rest/pug`)
- Secondary source: ChEMBL API (`https://www.ebi.ac.uk/chembl/api/data`)
- DB file: `db/fg_database.json`
- Always include `last_updated` timestamp in DB
- DB updater must be runnable as standalone: `python utils/db_updater.py`

### fg_database.json schema
```json
{
  "last_updated": "YYYY-MM-DD",
  "sources": ["pubchem", "chembl"],
  "functional_groups": {
    "Carboxylic acid": {
      "smarts": "C(=O)[OH]",
      "pubchem_cid": 280,
      "chembl_id": null,
      "description": "...",
      "known_target_classes": ["protease", "transporter"]
    }
  }
}
```

---

## Input formats (via OpenBabel)

Currently supported:
- CSV: col1 = compound name, col2 = SMILES

Planned (extend io_handler.py):
- SDF files
- MOL2 files
- InChI strings

Use OpenBabel for all format conversions. Do not write custom parsers.

---

## Output formats

- `fg_abundance_table.csv` — rows=functional groups, cols=compounds, values=integer counts
- `output/images/{compound_name}_fg.svg` — molecule with highlighted functional groups
- Future: PDF/HTML report for non-technical clients

---

## Key dependencies

```
rdkit          # Core cheminformatics
openbabel      # Format conversion
pandas         # Tables
requests       # API calls to PubChem / ChEMBL
matplotlib     # Plotting (future)
```

conda environment name: `chem_target`

---

## What NOT to do

- Do not use AutoDock or Glide — developer uses iGEMDOCK + SiMMap
- Do not merge constants into logic files
- Do not write to db/ manually — always via db_updater.py
- Do not over-engineer early — get it working first, then clean up
- Do not ignore OpenBabel for format handling — it is already installed

---

## Current status

| Component | Status |
|---|---|
| conda env + rdkit + openbabel | ✅ Done |
| v1 prototype (basic FG table) | ✅ Done |
| v2 prototype (readable names + SVG) | ✅ Done (monolithic) |
| Layered refactor | 🔲 Next |
| db_updater.py (PubChem + ChEMBL) | 🔲 Next |
| io_handler.py (CSV + SDF) | 🔲 Next |
| Target prediction layer | 🔲 Future |
| Report generator | 🔲 Future |

---

## Immediate tasks (in order)

1. Refactor v2 into layered architecture (constants / utils / main)
2. Build `utils/db_updater.py` — fetch from PubChem, fallback ChEMBL
3. Build `utils/io_handler.py` — CSV done, add SDF via OpenBabel
4. Review `data/*_pubchem_raw.json` output from explore_pubchem.py
   and finalize fg_database.json schema with developer
