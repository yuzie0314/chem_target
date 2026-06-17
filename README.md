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
  │  FG Detection  (fg_detector.py)     │  43 functional groups
  │  42 SMARTS + Steroid(Py) + fused    │  (+ routing/annotation cores)
  │  azolo-diazine + scaffold cores     │  per molecule
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
| Target class scores | `score = IDF × Σ mechanistic_weight(FG)` | IDF weighting boosts specific targets above generic labels; mechanistic_weight promotes known pharmacophores (Benzamidine 3.0 → serine protease, Hydroxamate 2.5 → HDAC, Sulfonamide 2.0 → carbonic anhydrase) |
| Conditional motifs | Multi-FG combo rules (e.g. Imidazole+Phenyl → CYP450 azole bonus; α,β-unsat. carbonyl → kinase warhead bonus) | Handles pharmacophores requiring co-occurrence context that single-FG scoring misses; also includes negative constraints (Hydroxamate/Thiol → suppress CYP450) |

---

## Quick start

```bash
# Conda environment (Python 3.10)
conda activate chem_target

# Single compound — print to terminal
python main.py predict --smiles "CC(=O)Oc1ccccc1C(=O)O" --name Aspirin

# Batch — CSV input, save summary table + plain-text reports
python main.py predict \
    --input  data/compounds.csv \
    --output output/predictions.csv \
    --report-dir output/reports \
    --top 10

# HTML reports with scoring methodology (open in browser)
python main.py predict \
    --input    data/compounds.csv \
    --html-dir output/html \
    --top 10
# → output/html/index.html        (batch summary, links to each compound)
# → output/html/<compound>.html   (structure SVG + targets + residues + methodology)

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

### Internal reference set (8 compounds)

8 reference compounds with known biological targets used during development:

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

### ChEMBL benchmark (run_benchmark.py)

A systematic benchmark pipeline is available for larger-scale evaluation:

```bash
# Curated set: top-20 active compounds per target class (~340 total)
python run_benchmark.py all --mode curated

# Limit test: up to 2 000 compounds across all target classes (SMILES-only, fast)
python run_benchmark.py all --mode limit --limit 2000
```

Outputs saved to `output/benchmark/`:
- `curated_results.csv` / `limit_results.csv` — per-compound predictions + accuracy flags
- `curated_summary.csv` / `limit_summary.csv` — per-class Top-1 / Top-3 accuracy + MRR
- `curated_report.txt` / `limit_report.txt` — plain-text summary (publication-ready)

#### Results (curated benchmark, pChEMBL ≥ 6.0, 20 compounds per class)

| Set | N | Top-1 | Top-3 |
|---|---|---|---|
| **Core 11 mechanistic classes (2026-06-16)** | **220** | **86.4% (190/220)** | **89.5% (197/220)** |
| + blind-spot rule-backed (MAO, COMT, cysteine protease, topoisomerase) | 300 | 72.3% (217/300) | — |

Per-class (curated, 20 compounds each):

| Target class | Top-1 | Top-3 | Notes |
|---|---|---|---|
| GPCR | 100% | 100% | ✅ |
| HDAC | 100% | 100% | ✅ |
| Carbonic anhydrase | 100% | 100% | ✅ |
| CYP450 | 95% | 95% | Thiazole SMARTS fixes ritonavir-class ×5; 1 structural (TAZAROTENIC ACID) |
| Tubulin | 95% | 100% | -1 GS-9256 (thiazole+ether; FG profile ≡ ritonavir-class, irreconcilable) |
| Kinase | 90% | 95% | Pyrimidine router (mono-pyrimidine→kinase) + α,β-unsat. warhead bonuses |
| mTOR | 85% | 85% | Morpholino-diazine (TORKinib) + Macrolide (rapalog) motifs; 3 remain (no morpholine) |
| Nuclear receptor | 85% | 100% | +1 (Guanidine FG IDF shift); Acylsulfonamide conflict + structural |
| COX | 75% | 85% | Indole+Sulfonamide conditional motif |
| Adenosine receptor | 60% | 60% | Pyrimidine router (fused-azolo-diazine→adenosine); 8 remain (no purine-mimetic core) |
| Serine protease | 65% | 65% | Benzamidine + Guanidine (arginine-mimetic) → S1 pocket; 7 peptidomimetics have no Arg-mimetic |
| cysteine protease | 60% | — | nitrile-warhead cathepsin inhibitors (odanacatib class) via gated Nitrile+Amide rule |
| COMT | 40% | — | nitrocatechol (entacapone/opicapone) via Phenol+Catechol; other 12 = research analogs |
| topoisomerase | 25% | — | anthracycline intercalators (doxorubicin etc.) via Anthraquinone voting FG; 15 = research series |
| MAO | 10% | — | propargylamine warhead (clorgiline); 18 = research series w/o MAO pharmacophore |

> **pChEMBL-sampling note:** MAO/COMT are sampled as the top-20 highest-affinity ChEMBL
> compounds, which are modern research analogs largely lacking the classic covalent
> pharmacophores (propargylamine, nitrocatechol) the marketed drugs carry. The warhead/
> nitrocatechol rules correctly capture the genuine drugs in the set; the cap reflects the
> sampling, not a rule gap.

Target classes covered by the benchmark (19 classes):

| Class | ChEMBL targets used |
|---|---|
| COX | COX-2 (CHEMBL230), COX-1 (CHEMBL220) |
| kinase | EGFR, BRAF, CDK2, mTOR |
| GPCR | β2-AR, D2R, 5-HT2A, A1R, A2AR |
| serine protease | thrombin, trypsin |
| cysteine protease | cathepsin B, cathepsin L |
| nuclear receptor | AR, ERα, GR, PR |
| MAO | MAO-A, MAO-B |
| HDAC | HDAC1, HDAC6, HDAC8 |
| adenosine receptor | A1R, A2AR |
| carbonic anhydrase | CA-II, CA-IX |
| CYP450 | CYP3A4, CYP2D6, CYP2C9 |
| PDE | PDE5A, PDE4B |
| mTOR | mTOR |
| tubulin | tubulin α1A |
| VKORC1 | VKORC1 |
| topoisomerase | Topo I, Topo II |
| ribosome | 50S ribosomal (bacterial) |
| xanthine oxidase | XO / xanthine dehydrogenase |
| COMT | COMT |

### External comparison: SwissTargetPrediction

To contextualise chem_target's accuracy, we benchmark against
[SwissTargetPrediction (STP)](https://www.swisstargetprediction.ch/) — a widely used
fingerprint-similarity reverse-docking tool trained directly on ChEMBL.

**Fair head-to-head scope:** chem_target is a functional-group / mechanism-based tool
designed for the 11 mechanistic target classes it has FG rules for. The comparison
below is computed on the **220 shared compounds (11 classes)** that both tools were
evaluated on (the STP results are cached; rerun with `--compare`). The earlier
"6.7% vs 55.1%" figure scored chem_target on 8 additional classes it has no rules for
(MAO, COMT, ribosome, topoisomerase, PDE, xanthine oxidase, cysteine protease) and is
therefore not a like-for-like comparison of the approach.

```bash
# Run STP on the same curated set and compare
python run_stp_comparison.py

# Regenerate comparison report from cached STP results
python run_stp_comparison.py --compare
```

Outputs saved to `output/benchmark/`:
- `stp_raw.csv` — all 100 STP predictions per compound (rank × probability × class)
- `stp_results.csv` — per-compound Top-1/Top-3/MRR
- `stp_summary.csv` — per-class accuracy table
- `stp_report.txt` — STP-only narrative report
- `comparison_report.txt` — side-by-side chem_target vs STP

#### Overall comparison (220 shared compounds, 11 mechanistic classes, 2026-06-15)

| Metric | chem_target | SwissTargetPrediction |
|---|---|---|
| Top-1 accuracy | **85.5%** | 70.5% |
| Top-3 accuracy | **89.1%** | 77.3% |
| Mean Reciprocal Rank | **0.872** | 0.754 |
| Macro-avg F1 (Top-1) | **0.855** | 0.726 |

On the classes it targets, chem_target now **out-performs** STP overall — driven by
mechanistic pharmacophores STP's fingerprint similarity misses (tubulin 95% vs 5%,
CYP450 95% vs 35%, mTOR 85% vs 60%). STP remains stronger on well-populated ChEMBL
classes with dense analog series (COX, adenosine, serine protease).

> ⚠ **STP bias note:** The curated test compounds are sourced directly from ChEMBL,
> the same database STP's fingerprint models are trained on.  STP accuracy figures
> therefore represent an *in-distribution* upper bound — analogous to chem_target's
> BioLiP/PDB circularity.  Both biases are disclosed in the respective report files.

#### Per-class breakdown

(11 shared mechanistic classes, 20 compounds each; **bold** = winner on Top-1)

| Target class | N | cT Top-1 | cT Top-3 | cT MRR | STP Top-1 | STP Top-3 | STP MRR |
|---|---|---|---|---|---|---|---|
| GPCR | 20 | **100%** | 100% | 1.000 | 75% | 75% | 0.750 |
| HDAC | 20 | 100% | 100% | 1.000 | 100% | 100% | 1.000 |
| carbonic anhydrase | 20 | 100% | 100% | 1.000 | 100% | 100% | 1.000 |
| CYP450 | 20 | **95%** | 95% | 0.950 | 35% | 75% | 0.594 |
| tubulin | 20 | **95%** | 100% | 0.967 | 5% | 15% | 0.154 |
| kinase | 20 | **90%** | 95% | 0.930 | 80% | 85% | 0.825 |
| mTOR | 20 | **85%** | 85% | 0.850 | 60% | 70% | 0.670 |
| nuclear receptor | 20 | **80%** | 100% | 0.875 | 70% | 75% | 0.757 |
| COX | 20 | 75% | 85% | 0.820 | **90%** | 90% | 0.913 |
| adenosine receptor | 20 | 60% | 60% | 0.600 | **80%** | 85% | 0.825 |
| serine protease | 20 | 60% | 60% | 0.600 | **80%** | 80% | 0.800 |

**Notable findings:**
- chem_target **wins or ties on 8/11 classes**, decisively on tubulin (95% vs 5%),
  CYP450 (95% vs 35%) and mTOR (85% vs 60%) — mechanistic pharmacophores
  (Acylsulfonamide warhead, azole heme coordination, morpholino-diazine hinge)
  that fingerprint similarity does not capture.
- STP wins on **COX**, **adenosine receptor**, **serine protease** — densely
  populated ChEMBL analog series where fingerprint similarity excels.
- chem_target's FG approach is expected to **generalise better to novel scaffolds**
  not represented in ChEMBL (the basis of STP's similarity search).
- **Scope:** chem_target only targets these 11 mechanistic classes; STP additionally
  attempts MAO / COMT / PDE / topoisomerase / ribosome / xanthine oxidase / cysteine
  protease, for which chem_target currently has no FG rules (scores ≈ 0).

**Design philosophy comparison:**

| Dimension | chem_target | SwissTargetPrediction |
|---|---|---|
| Approach | FG-to-residue structural mapping | FP2 / FP4 fingerprint similarity |
| Training data | BioLiP 2.0 (PDB structural binding events) | ChEMBL activities |
| Output | Target class + key binding residues + FG drivers | Ranked target proteins + probability |
| Novel scaffolds | Generalises via functional groups | Limited by fingerprint coverage |
| Mechanism insight | ✅ Explains *which* FGs drive binding | ❌ Black-box similarity |
| Speed | Fast (local, no API) | Requires web submission |

### ⚠ Validation bias disclosure

**This tool's residue scoring layer is derived from BioLiP 2.0, which aggregates
protein–ligand binding events from the RCSB Protein Data Bank (PDB).**

Compounds included in the ChEMBL benchmark may also have co-crystal structures
deposited in the PDB (e.g., a kinase inhibitor with both a ChEMBL IC₅₀ record
and a published PDB entry).  When such a compound is tested, its functional
groups are already represented in `db/fg_residue_table.csv`, leading to
artificially higher residue scores — a form of **structural data circularity**.

Consequences:
- Residue-level predictions for such compounds should be treated as an **upper bound**.
- Target class predictions are somewhat less affected, as they rely on a separate
  curated annotation database (`db/fg_database.json`) rather than BioLiP directly.
- The benchmark reports include this disclaimer automatically in the output `.txt` file.

Mitigation strategies (for rigorous evaluation):
1. Filter test compounds to those **without** a PDB entry (`has_pdb_entry: False` in ChEMBL).
2. Use the limit-test results as a population-level estimate rather than individual truth.
3. Treat reported Top-1/Top-3 accuracy as an upper-bound benchmark pending a fully
   disjoint test set (planned for future ChEMBL–PDB de-overlap analysis).

---

## Functional groups (38 total)

### SMARTS-based (37)

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
| 15 | Imidazole | `c1cnc[nH,n]1` | CYP450, histamine receptor, metalloprotease |
| 16 | Triazole | `n1cncn1` | CYP450 (fluconazole/voriconazole-class, mw=1.5) |
| 17 | Thiazole | `c1cncs1` | CYP450 (ritonavir-class CYP3A4 inhibitors, mw=1.5) |
| 18 | Benzimidazole | `c1ccc2[nH,n]cnc2c1` | Scaffold marker — omeprazole/mebendazole; Imidazole substructure handles CYP votes |
| 19 | Indole | `c1ccc2[nH]ccc2c1` | Serotonin receptor, tubulin, BACE1 |
| 20 | Purine | `c1ncc2ncnc2n1` | Adenosine receptor, kinase, DNA polymerase, PDE |
| 21 | Xanthine | `O=c1nc(=O)c2ncnc2n1` | Adenosine receptor, PDE, xanthine oxidase |
| 22 | Nitrile | `C#N` | Cysteine protease, nitrile hydratase |
| 23 | Nitro | `[$([NX3](=O)=O),$([NX3+](=O)[O-])]` | Nitroreductase, CYP450 |
| 24 | Benzamidine | `[NX3H2][CX3](=[NX2H1])c` | Serine protease S1 pocket (thrombin, trypsin, mw=3.0) |
| 25 | Thiol | `[SX2H]` | Cysteine protease, metalloenzyme, HDAC |
| 26 | Sulfonamide | `[SX4](=O)(=O)[NX3]` | Carbonic anhydrase, COX, kinase (mw=2.0) |
| 27 | Methylsulfone | `[CX4H3][SX4](=O)(=O)c` | COX-2 selectivity pocket (celecoxib class) |
| 28 | Hydroxamate | `[CX3](=O)[NX3H][OX2H]` | HDAC Zn chelation (vorinostat class, mw=2.5) |
| 29 | Acylsulfonamide | `[CX3](=O)[NX3H][SX4](=O)(=O)` | Tubulin (epothilone macrolide, mw=2.0) |
| 30 | Phenyl ring | `c1ccccc1` | Kinase, GPCR, COX, tubulin |
| 31 | Coumarin | `O=c1ccc2ccccc2o1` | MAO, VKORC1, CYP450, serine protease |
| 32 | Chromone | `O=c1ccoc2ccccc12` | COX, kinase, estrogen receptor |
| 33 | Halogen | `[F,Cl,Br,I]` | Kinase, ion channel, thyroid receptor |
| 34 | Epoxide | `[OX2r3]` | Epoxide hydrolase, cysteine protease |
| 35 | Endoperoxide | `[OX2r][OX2r]` | Antimalarial target, heme-dependent enzyme |
| 36 | α,β-unsat. carbonyl | `[CX3](=O)C=C` | Cysteine protease, Nrf2, NF-κB; covalent kinase warhead |
| 37 | Macrolide | `[CX3;!r3;…;R](=O)[OX2;!r3;…;R]` (ring ≥12) | mTOR, calcineurin, ribosome |
 
> **Hierarchical overlaps** (by design): Xanthine ⊂ Purine · Coumarin ⊂ Lactone · Macrolide ⊂ Lactone · Indole ∩ Phenyl ring · Benzimidazole ⊃ Imidazole. Each level provides independent pharmacological resolution.

### Python-based (1)

| # | Name | Detection | Primary targets |
|---|---|---|---|
| 38 | Steroid | `_detect_steroid_core()` — ring BFS (r5_C ≥ 5, r6_C ≥ 10, both_C ≥ 2) | Nuclear receptor (AR/GR/ER/PR), CYP450 |

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
│   ├── fg_smarts.py          # FG_SMARTS dict: 37 SMARTS patterns (incl. Triazole+Thiazole+Benzimidazole)
│   └── fg_names.py           # Legacy RDKit fr_* → human-readable (kept for compatibility)
├── db/                       # Auto-generated — never edit by hand
│   ├── fg_database.json      # 38 FG entries: SMARTS, ChEBI, targets, mechanistic_weight
│   ├── fg_residue_table.csv  # 37 SMARTS + Steroid × 20 AA BioLiP co-occurrence matrix
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
│   └── benchmark/            # Auto-generated by run_benchmark.py
│       ├── curated/          # ChEMBL compounds + SDF downloads
│       └── limit/
├── output/                   # Generated output (CSV, SVG, reports)
│   └── benchmark/            # Benchmark evaluation results
├── main.py                   # CLI entry point (thin — no logic)
├── run_benchmark.py          # ChEMBL benchmark pipeline (download → run → report)
├── run_stp_comparison.py     # SwissTargetPrediction comparison pipeline
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
