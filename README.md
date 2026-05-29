# chem\_target

**Compound-to-target prediction tool** (reverse docking / target fishing).

Input a SMILES string or a CSV/SDF file of molecules.  
Output: predicted enzyme/protein target classes and key binding residues,
ranked by a two-stage weighted score derived from structural binding data
and pharmacological annotations.

---

## Pipeline overview

```
Input (SMILES / CSV / SDF)
        │
        ▼
  ┌─────────────────────────────────────┐
  │  FG Detection  (fg_detector.py)     │  32 functional groups
  │  31 SMARTS patterns + Steroid(Py)   │  per molecule
  └─────────────────────────────────────┘
        │
        ▼
  ┌─────────────────────────────────────┐   db/fg_residue_table.csv
  │  Residue Scoring                    │ ◄─ BioLiP 2.0 (46 114 events,
  │  predict_residues()                 │    6 002 unique ligands)
  │  z-score normalised sum             │
  └─────────────────────────────────────┘
        │
        ▼
  ┌─────────────────────────────────────┐   db/fg_database.json
  │  Target Class Voting                │ ◄─ ChEBI / PubChem / manual
  │  predict_target_classes()           │    annotations per FG
  │  IDF-weighted votes                 │
  └─────────────────────────────────────┘
        │
        ▼
  Text report / Summary CSV / SVG images
```

### Scoring design

| Layer | Method | Why |
|---|---|---|
| Residue scores | Per-FG z-score normalisation → sum | Prevents high-frequency generic FGs (e.g. Hydroxyl: ~11 000 GLY events) from drowning out specific pharmacophores (e.g. Steroid: ~1 400 PHE events) |
| Target class scores | `score = votes × log(N_FGs / N_FGs_with_this_target)` | IDF weighting boosts specific targets (VKORC1, tubulin, antimalarial target) above generic labels (kinase, GPCR) even with fewer raw votes |

---

## Quick start

```bash
# Conda environment (Python 3.10)
conda activate chem_target

# Single compound — print to terminal
python main.py predict --smiles "CC(=O)Oc1ccccc1C(=O)O" --name Aspirin

# Batch — CSV input, save summary table + individual reports
python main.py predict \
    --input  data/compounds.csv \
    --output output/predictions.csv \
    --report-dir output/reports \
    --top 10

# Functional group abundance table + SVG images
python main.py fg \
    --input  data/compounds.csv \
    --output fg_abundance_table.csv \
    --images output/images
```

### Input CSV format

```
Name,SMILES
Aspirin,CC(=O)Oc1ccccc1C(=O)O
Caffeine,Cn1cnc2c1c(=O)n(C)c(=O)n2C
Warfarin,CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O
```

### Output summary CSV columns

| Column | Description |
|---|---|
| `compound` | Display name |
| `smiles` | Input SMILES |
| `n_fgs` | Number of FGs detected |
| `fgs` | FG names (comma-separated) |
| `target_1/2/3` | Top-3 predicted target classes |
| `score_1/2/3` | IDF-weighted scores |
| `votes_1/2/3` | Raw vote counts |
| `residue_1/2/3` | Top-3 binding residues (3-letter code) |
| `res_score_1/2/3` | Z-score-normalised residue scores |
| `warning` | Error message if SMILES invalid or no FGs found |

---

## Validation results

8 reference compounds with known biological targets:

| Compound | Known target | Top prediction | Score |
|---|---|---|---|
| Aspirin | COX-1/2 | COX | 4.93 (3 votes) |
| Caffeine | Adenosine receptor / PDE | Adenosine receptor | 5.48 (2 votes) |
| Warfarin | VKORC1 / CYP450 | VKORC1 + anticoagulant target | 3.43 |
| Quercetin | COX / kinase / ER | Estrogen receptor | 6.14 (3 votes) |
| Artemisinin | Antimalarial / heme | Antimalarial target | 3.43 |
| Erythromycin | Ribosome / mTOR | Ribosome + mTOR | 3.43 each |
| Vincristine | Tubulin | Tubulin | 5.48 (2 votes) |
| Testosterone | Nuclear receptor (AR/GR) | Androgen / glucocorticoid receptor | 3.43 each |

---

## Functional groups (32 total)

### SMARTS-based (31)

| # | Name | SMARTS | Primary targets |
|---|---|---|---|
| 1 | Carboxylic acid | `C(=O)[OH]` | Protease, COX, ACE inhibitor |
| 2 | Ester | `[#6][CX3;!R](=O)[OX2H0;!R][#6]` | Esterase, lipase, COX |
| 3 | Amide | `[CX3](=O)[NX3;H1,H2]` | Protease, kinase |
| 4 | Lactone | `[#6X3;R](=O)[#8X2;R]` | HDAC, topoisomerase I, esterase |
| 5 | Ketone | `[#6][#6X3H0](=O)[#6]` | Oxidoreductase, 17β-HSD |
| 6 | Aldehyde | `[CX3H1](=O)` | ALDH, ADH, protease |
| 7 | Hydroxyl | `[OX2H;!$(Oc);!$(OC=O)]` | Kinase, phosphatase, transporter |
| 8 | Phenol | `[OX2H]c` | COX, COMT, estrogen receptor |
| 9 | Catechol | `c1cc([OX2H])c([OX2H])cc1` | COMT, MAO, dopamine receptor |
| 10 | Ether | `[OX2;!$(O=*);!$([OX2H]);!$(OC=O)]([#6])[#6]` | CYP450 |
| 11 | Methylenedioxy | `c1ccc2c(c1)OCO2` | CYP450, MAO |
| 12 | Primary amine | `[NX3H2;!$(NC=O);!$(NS=O);!$(Nc)]` | MAO, GPCR, transporter |
| 13 | Secondary amine | `[NX3H1;!$(NC=O);!$(NS=O);!$(Nc)]` | MAO, GPCR, ion channel |
| 14 | Tertiary amine | `[NX3H0;!$(N=*);!$(NC=O);!$(NS=O);!$(Nc)]` | GPCR, nicotinic receptor |
| 15 | Imidazole | `c1cnc[nH]1` | CYP450, histamine receptor, metalloprotease |
| 16 | Indole | `c1ccc2[nH]ccc2c1` | Serotonin receptor, tubulin, BACE1 |
| 17 | Purine | `c1ncc2ncnc2n1` | Adenosine receptor, kinase, DNA polymerase, PDE |
| 18 | Xanthine | `O=c1nc(=O)c2ncnc2n1` | Adenosine receptor, PDE, xanthine oxidase |
| 19 | Nitrile | `C#N` | Cysteine protease, nitrile hydratase |
| 20 | Nitro | `[$([NX3](=O)=O),$([NX3+](=O)[O-])]` | Nitroreductase, CYP450 |
| 21 | Thiol | `[SX2H]` | Cysteine protease, metalloenzyme, HDAC |
| 22 | Sulfonamide | `[SX4](=O)(=O)[NX3]` | Carbonic anhydrase, COX, kinase |
| 23 | Phenyl ring | `c1ccccc1` | Kinase, GPCR, COX, tubulin |
| 24 | Coumarin | `O=c1ccc2ccccc2o1` | MAO, VKORC1, CYP450, serine protease |
| 25 | Chromone | `O=c1ccoc2ccccc12` | COX, kinase, estrogen receptor |
| 26 | Halogen | `[F,Cl,Br,I]` | Kinase, ion channel, thyroid receptor |
| 27 | Epoxide | `[OX2r3]` | Epoxide hydrolase, cysteine protease |
| 28 | Endoperoxide | `[OX2r][OX2r]` | Antimalarial target, heme-dependent enzyme |
| 29 | α,β-unsat. carbonyl | `[CX3](=O)C=C` | Cysteine protease, Nrf2, NF-κB |
| 30 | Macrolide | `[CX3;!r3;…;R](=O)[OX2;!r3;…;R]` (ring ≥12) | mTOR, calcineurin, ribosome |
| 31 | Methylenedioxy | `c1ccc2c(c1)OCO2` | CYP450, MAO |

> **Hierarchical overlaps** (by design): Xanthine ⊂ Purine · Coumarin ⊂ Lactone · Macrolide ⊂ Lactone · Indole ∩ Phenyl ring. Each level provides independent pharmacological resolution.

### Python-based (1)

| # | Name | Detection | Primary targets |
|---|---|---|---|
| 32 | Steroid | `_detect_steroid_core()` — ring BFS (r5_C ≥ 5, r6_C ≥ 10, both_C ≥ 2) | Nuclear receptor (AR/GR/ER/PR), CYP450 |

> **Why Python?** RDKit's `rN` SMARTS primitive uses the Smallest Set of Smallest Rings (SSSR). In the 6-6-6-5 steroid tetracycle, C-D ring junction atoms are assigned only to the smallest ring (r5), making `[r5;r6]` always fail. `IsAtomInRingOfSize()` is not SSSR-dependent and correctly identifies junction atoms.

---

## Data sources

| Source | Role | Size |
|---|---|---|
| **BioLiP 2.0** (`db/BioLiP_nr.txt.gz`) | FG × residue co-occurrence table | 46 114 binding events, 6 002 unique ligands |
| **RCSB CCD** (via REST API) | SMILES for BioLiP ligand 3-letter codes | 6 002 entries cached in `db/ccd_smiles_cache.json` |
| **PubChem / ChEBI** | FG metadata (CID, ChEBI ID, description) | `db/fg_database.json` |
| **RCSB PDB** (via REST API) | 3D residue–ligand poses | 3 222 records in `db/residue_3d_poses.json`, 86 SDF files |

---

## Project structure

```
chem_target/
├── constants/
│   ├── fg_smarts.py          # FG_SMARTS dict: 31 SMARTS patterns
│   └── fg_names.py           # Legacy RDKit fr_* → human-readable (kept for compatibility)
├── db/                       # Auto-generated — never edit by hand
│   ├── fg_database.json      # 32 FG entries: SMARTS, ChEBI, targets, descriptions
│   ├── fg_residue_table.csv  # 31 FG × 20 AA BioLiP co-occurrence matrix
│   ├── ccd_smiles_cache.json # 6 002 RCSB CCD SMILES entries
│   ├── residue_3d_poses.json # 3 222 Cα + ligand centroid + distance records
│   ├── local_env/*.sdf       # 86 representative FG–residue complex SDFs
│   ├── pharmacophore_stats.json
│   └── BioLiP_nr.txt.gz      # Raw BioLiP source (not tracked in git)
├── utils/
│   ├── fg_detector.py        # detect_smarts(), _detect_steroid_core(), _PYTHON_DETECTORS
│   ├── target_predictor.py   # predict(), predict_residues(), predict_target_classes()
│   ├── interaction_analyzer.py  # BioLiP → fg_residue_table.csv builder
│   ├── pose_extractor.py     # 3D pose → residue_3d_poses.json + SDF builder
│   ├── db_updater.py         # PubChem / ChEMBL → fg_database.json
│   ├── io_handler.py         # CSV / SDF reading
│   └── visualizer.py         # RDKit SVG output
├── data/                     # User input files
├── output/                   # Generated output (CSV, SVG, reports)
├── main.py                   # CLI entry point (thin — no logic)
└── README.md
```

---

## Rebuilding the database

Run these in order after changing `constants/fg_smarts.py` or `utils/fg_detector.py`:

```bash
# 1. Rebuild FG × residue co-occurrence table (~5 min, uses SMILES cache)
python utils/interaction_analyzer.py --local db/BioLiP_nr.txt.gz

# 2. Rebuild 3D poses, SDF files, pharmacophore stats (~20 min, downloads PDBs)
python -u utils/pose_extractor.py --top-per-fg 3
```

> **Windows / Unicode note:** use the conda environment Python directly:  
> `$env:PYTHONIOENCODING="utf-8"; & "C:\...\envs\chem_target\python.exe" -u utils/interaction_analyzer.py ...`

---

## Dependencies

```
rdkit          >= 2023      Core cheminformatics
openbabel      >= 3.1       SDF / MOL2 format conversion
pandas         >= 1.5       Tabular data
requests       >= 2.28      PubChem / ChEMBL / RCSB API calls
biopython      >= 1.79      PDB parsing (pose_extractor)
numpy          >= 1.23      Numerical operations
```

Install via conda:

```bash
conda create -n chem_target python=3.10
conda activate chem_target
conda install -c conda-forge rdkit openbabel biopython pandas requests numpy
```

---

## Architecture rules

- **`constants/`** — static lookup tables only; no logic
- **`db/`** — auto-generated only; never edit by hand
- **`main.py`** — argument parsing and dispatch only; no business logic
- **`utils/`** — pure functions where possible; no side effects outside designated builders

Adding a new Python-based FG detector (like Steroid):
1. Implement in `utils/fg_detector.py`, add to `_PYTHON_DETECTORS`
2. Add SMARTS-free entry to `db/fg_database.json` (`"smarts": null`)
3. Verify `utils/interaction_analyzer.py` and `utils/pose_extractor.py` import the detector
4. Rebuild both DB files (see above)

---

## Known limitations

| Limitation | Affected compounds | Notes |
|---|---|---|
| No β-lactam FG | Penicillin, cephalosporins | 4-membered lactam ring not detected |
| No quinoline / isoquinoline FG | Chloroquine, quinine | Bicyclic N-heterocycle coverage gap |
| No guanidine FG | Metformin, arginine analogs | Strongly basic, unique salt-bridge profile |
| BioLiP bias toward fragment-like ligands | Large natural products | Crystal structure ligands tend to be small; macrolide and steroid counts lower than in vivo relevance |
| Layer 2 custom DB not yet built | Activity-based targets | Planned: ChEMBL IC₅₀/Kᵢ aggregated by FG profile |

---

## Roadmap

- [ ] HTML / PDF report for non-technical clients
- [ ] Layer 2 custom interaction database (ChEMBL activity-based)
- [ ] Batch SDF input with property-based name extraction
- [ ] Additional scaffolds: β-lactam, quinoline, guanidine
- [ ] Web UI / API endpoint (SaaS phase)
- [ ] Integration with iGEMDOCK / SiMMap for structure-based second-stage filter

---

## Author

Yu-Jen Chang — pharmacist, computational chemist, data scientist.  
Primary docking tools: iGEMDOCK, SiMMap.  
Stack: Python 3.10 · RDKit · OpenBabel · pandas · conda.
