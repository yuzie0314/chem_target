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
│   ├── fallback_3d.py         # 3D-fallback (Fallback3D, ProLIFFallback: pH-protonate→per-reference smina dock→ProLIF IFP→Jaccard, build_override) — lazy heavy deps
│   ├── build_prolif_reference.py  # env check + ProLIF reference-IFP builder (check|build) → db/prolif_reference_ifp.json
│   ├── interaction_analyzer.py  # BioLiP → fg_residue_table.csv
│   ├── pose_extractor.py      # 3D pose extractor → residue_3d_poses.json
│   ├── db_updater.py          # PubChem/ChEMBL → fg_database.json
│   ├── io_handler.py          # CSV/SDF input parsing
│   ├── report_generator.py    # HTML individual + batch reports
│   ├── mlflow_tracker.py      # Opt-in MLflow logging (SQLite store) — lazy import, no-op by default
│   └── visualizer.py          # RDKit SVG output
├── data/           # User input files
├── output/         # Generated output — CSV, SVG, HTML, reports
├── run_benchmark.py   # 11-class × 20-compound benchmark pipeline (--mlflow opt-in)
├── backfill_mlflow.py # One-off: backfill pre-DVC milestones as reconstructed MLflow runs
├── eval_holdout.py    # Held-out generalisation eval: tuning vs zero-overlap held-out gap
├── tests/             # Stdlib unittest regression net (locks the headline; no pytest dep)
└── main.py         # Entry point only. Thin. No logic here.
```

DVC-tracked db artifacts (not in git, see "Reproducibility" below): `db/BioLiP_nr.txt.gz`,
`db/local_env/`, `db/prolif_pdb/`, `db/residue_3d_poses.json`, `db/ccd_smiles_cache.json`.

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

dvc            # Data version control for large db/ artifacts

mlflow         # Benchmark experiment tracking (opt-in)

```



conda environment name: `chem\_target`



\---



\## Reproducibility: DVC + MLflow



**DVC (local cache, no remote yet)** — 5 large/binary `db/` artifacts are
DVC-tracked, not in git (git holds only the `.dvc` pointer files):
`BioLiP_nr.txt.gz` (14M), `local_env/` (7.9M), `prolif_pdb/` (5.9M),
`residue_3d_poses.json` (1.9M), `ccd_smiles_cache.json` (420K).
Small out-of-box artifacts stay in git (`fg_database.json`,
`fg_residue_table.csv`, `target_class_aliases.json`, `prolif_reference_ifp.json`).
`dvc status` to check; `dvc checkout` to restore working copies from cache.
**No off-machine remote yet** → local cache is the only copy of the
non-regenerable data; add a `dvc remote` for backup before relying on it.
**No `dvc.yaml` pipeline** (regeneration scripts are network-bound /
non-deterministic — deliberately skipped per "don't over-engineer").

**MLflow (opt-in, local SQLite store)** — `utils/mlflow_tracker.py` logs to
`mlflow.db` + `mlartifacts/` (both gitignored; MLflow 3 retired the file store).
- `python run_benchmark.py report --mlflow` (or `all --mlflow`) records the run:
  full-16-class + **core-11-subset** Top-1/Top-3, per-class Top-1, MRR,
  macro/weighted F1, git SHA/dirty, + report & summary as artifacts.
  Default (no flag) is byte-identical — mlflow is never imported unless asked.
- `core11_*` metrics = the headline core-11 subset basis (190/220 when the
  backfilled runs were recorded; currently 185/220 after the 2026-07-03
  de-overfitting pass); plain `top1_acc` on the live run is the wider
  320/16-class basis (do NOT compare the two series directly).
- `python backfill_mlflow.py` = one-off; 17 pre-DVC milestones (59.1%→86.4%)
  reconstructed from commit messages, tagged `reconstructed=true`. **Run once
  only** (re-running duplicates; clear `mlflow.db` first). Pre-DVC dataset
  snapshots are GONE (gitignored + overwritten), so older points are
  reconstructed, not rerunnable — the first fully reproducible run is HEAD.
- View: `mlflow ui --backend-store-uri sqlite:///mlflow.db`.



\---



\## What NOT to do

\- Do not use AutoDock or Glide — developer uses iGEMDOCK + SiMMap
\- Do not merge constants into logic files
\- Do not write to db/ manually — always via interaction\_analyzer.py / pose\_extractor.py / db\_updater.py
\- Do not over-engineer early — get it working first, then clean up
\- Do not ignore OpenBabel for format handling — it is already installed



\---



\## Current status (headline)

- **Core 11-class Top-1 (TUNING set): 185/220 = 84.1%** (mechanistic classes) — branch `dev/validation`.
  - ⚠️ **This is a resubstitution number.** All conditional rules were tuned against this exact curated set (now frozen in git). It is an **upper bound**, not a generalisation estimate.
  - **PURE HELD-OUT Top-1: 585/936 = 62.5%** (`data/benchmark/holdout/compounds.csv` = `limit` minus curated by id AND canonical SMILES, zero overlap; leak-guarded in eval_holdout.py) → **+21.6-pt overfitting gap**. Hand-tuned classes collapse out-of-sample (CYP450 70→27%, serine protease 65→29%, COX 75→38%, tubulin 95→41%); pharmacophore-specific classes hold (CA 100→95%, GPCR 100→94%). Reproduce: `python eval_holdout.py`.
  - **2026-07-03 de-overfitting pass** (ablation-driven, see "De-overfitting" below): removed the ID-tuned CYP450 aryl-halide rules (−5 tuning CYP450, +0 held-out — pure memorisation) and re-keyed the cathepsin rule on a *non-aryl* nitrile warhead (+5 held-out, mostly NR; all 12 curated cathepsin TPs kept). Gap 24.4→21.6pt. Was 190/220=86.4%.
- Blind-spot rule-backed (extended set): MAO 2/20, COMT 8/20, cysteine protease 12/20, topoisomerase 5/20.
- The rule layer is complete; ProLIF 3D-fallback runs end-to-end but is NOT auto-registered (zero regression) and recovers 0/7 SP misses (see the 3D-fallback section).
- Per-class detail → `doc/benchmark_results.md`; the "what's built" checklist lives in `git log` (the rules/decisions that matter are in the sections that follow — mechanistic_weight, conditional rules, negative results, structural limits).

---

\## Benchmark per-class results

Full per-class table + Notes (the per-class reasoning / structural-limit provenance), the
13-class extended-set pChEMBL-bias note, and the run_stp_comparison data-bug fix → **`doc/benchmark_results.md`**.
The AI-critical "do NOT" facts are also distilled in "Known structural limitations" + "Confirmed
negative results" below.

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

`Non-aryl nitrile` — SMARTS `[!c;!#1][CX2]#[NX1]` in `_WARHEAD_ANNOTATIONS` (routing-only, not in fg_database → no IDF shift). A nitrile NOT bonded to an aromatic ring = the electrophilic cathepsin thioimidate warhead (aliphatic/heteroatom-bound), as opposed to an inert aryl nitrile. Consumed only by the cysteine-protease rule (replaced the generic `Nitrile` FG there, 2026-07-03).

---

\## Conditional scoring rules (utils/target_predictor.py)

All rules are pre-IDF bonuses (multiplied by IDF before adding to final score).

| Rule | Condition | Target | Bonus | Rationale |
|---|---|---|---|---|
| CYP450 azole | **free heme azole** (`_has_free_heme_azole`) + {Phenyl/Ether/Halogen}, Ketone only if no Amide/TertAmine, no Purine/αβunsat/Sulfonamide | cytochrome P450 | +2.0 | Azole heme-Fe coordination (fluconazole/voriconazole/ketoconazole/ritonavir class). Free azole = Triazole OR Thiazole OR (Imidazole & not Benzimidazole), AND not Fused-azolo-diazine — excludes purine-mimetic fused cores (adenosine/kinase) while keeping voriconazole (free triazole + separate fluoropyrimidine). **The ONLY CYP450 rule** — the aryl-halide/amide-halide/ether-amine COOH sub-rules were removed 2026-07-03 (ID-tuned memorisation: −5 tuning, +0 held-out) |
| COX indole-sulfonamide | Indole + Sulfonamide | COX | +2.0 | Indole scaffold + COX-2 selectivity pocket |
| mTOR macrolide | Macrolide, no Thiol/αβunsat/Acylsulfonamide | mTOR | +2.0 | Rapamycin-class allosteric FKBP12 binding |
| Adenosine Purine | Purine present | adenosine receptor | +0.5 | Purine is the defining adenosine scaffold |
| Kinase αβunsat warhead | α,β-unsat. carbonyl present | kinase | +0.5 | Covalent Michael acceptor warhead (EGFR) |
| Kinase sulfonamide-amine | Sulfonamide + TertAmine | kinase | +2.0 | Kinase linker hijacked by CA (Sulfonamide mw=2.0) |
| MAO warhead | (Propargylamine OR Hydrazine), no Sulfonamide/Nitrile/αβunsat | MAO | +2.5 | Irreversible MAO inhibitor: propargylamine→FAD adduct (selegiline/clorgiline) or hydrazine (phenelzine). Exclusions prevent CA/covalent-kinase false positives. Markers are routing-only (`_WARHEAD_ANNOTATIONS`, not in fg_database → no IDF shift) |
| Cysteine protease nitrile | **Non-aryl nitrile** + Amide, no αβunsat/Pyrimidine/Quinazoline/Pyrrolopyrimidine/Fused-azolo/Sulfonamide | cysteine protease | +2.5 | Peptidomimetic nitrile warhead → thioimidate with catalytic Cys25 (odanacatib/cathepsin-K class). **Keys on the "Non-aryl nitrile" routing marker (2026-07-03), not the generic Nitrile FG** — an aryl nitrile (enobosarm/aryl-CN androgen ligands) is inert, not a warhead; the switch stops the rule mis-grabbing aryl-nitrile NR/HDAC compounds (held-out +5) while keeping all 12 curated cathepsin TPs. Exclusions strip kinase-hinge / covalent-kinase / CA contexts |

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
3. **CYP450 70% (14/20)** (was 95% until 2026-07-03): the ID-tuned aryl-halide/amide-halide/ether-amine COOH rules were removed (memorised specific ChEMBL ids; +0 on held-out — see "De-overfitting" in the headline). Only the mechanistic azole heme-Fe rule remains (Thiazole SMARTS still fixes ritonavir-class ×5). The 5 lost compounds are non-azole lipophilic substrates with no generalisable CYP pharmacophore — **do NOT re-add halogen/COOH-combo rules to chase them** (that is exactly the memorisation that was removed).
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

\### 3D-fallback layer — runs end-to-end, but DO NOT register for a benchmark gain

Full design / builder / tune-list → **`doc/3d_fallback.md`**; user-facing write-up → README
"3D interaction-fingerprint fallback". The load-bearing rules:

- **Confidence gate** (`target_predictor.assess_confidence` + `register_fallback_3d` +
  `_finalize_with_fallback`): low/none-confidence predictions route to a pluggable 3D hook; default
  `_stub_fallback_3d` → None → **zero regression**. The gate is PARTIAL — 63/105 misses are
  high-confidence-wrong and unreachable by it (needs a parallel path + meta-reconciler, regression risk).
- **ProLIF fallback** (`utils/fallback_3d.py`) is fully implemented and runs end-to-end (smina +
  PDBFixer + pH 7.4 protonation + per-reference docking + chain-stripped IFP keys; 9-co-crystal SP
  reference library). **It is NOT auto-registered** — predict/benchmark are byte-identical (core
  190/220 unchanged).
- **⚠ EMPIRICAL: recovers 0/7 SP misses → NO accuracy gain (2026-06-18).** Best IFP Jaccard for the real
  misses is 0.33–0.60 (below the 0.6 fire threshold; flipping top-1 needs ≈0.72–0.82). Positive controls
  (benzamidine/rivaroxaban →1.0) only score high because they ARE in the reference set (self-docking).
  Confirms the 7 SP misses are real structural misses, not a pattern gap. **Do NOT register
  `ProLIFFallback` expecting a benchmark gain**; threshold-lowering won't help (fires on negatives).
  Payoff needs a far broader reference library / docking into the TRUE target receptor — see doc.

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
4. ~~ProLIF 3D-fallback~~ ✅ RUNS END-TO-END (2026-06-18): smina + PDBFixer protein-prep + pH 7.4
   query protonation + per-reference docking + chain-stripped IFP keys; 9-co-crystal SP reference
   library (4 peptidomimetic FXa added). Validated: benzamidine 1.0, rivaroxaban 1.0, ibuprofen 0.25.
   Still NOT auto-registered → zero regression. Next: per-ref docking cost, more SP co-crystals, then
   measure net gain of `register_fallback_3d(ProLIFFallback())` on the SP misses vs 41 low-conf-hit risk.

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

