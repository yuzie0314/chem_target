# 3D-fallback layer — full design reference

Detailed reference for the structure-based second-opinion layer. CLAUDE.md keeps only
the load-bearing rules (zero-regression invariant + the "do NOT register expecting a
benchmark gain" finding); this file has the mechanism, the builder, and the tune-list.
The user-facing write-up of the negative result is in `README.md` ("3D
interaction-fingerprint fallback (ProLIF)").

## Confidence gate (routing skeleton)

`utils/target_predictor.py`: `assess_confidence()` + `register_fallback_3d()` +
`_finalize_with_fallback()`. Predictions get a `confidence` tag (high/low/none) from
observable signal only (top1 score < 3.0 OR top1−top2 margin < 0.75 → low; empty votes →
none). Low/none route to a pluggable 3D hook; default `_stub_fallback_3d` returns None →
**zero regression** (verified: core stays unchanged). Plug real models via
`register_fallback_3d`.

**Gate-selectivity finding (320-cpd benchmark)** — the confidence gate is SAFE but PARTIAL:
- Routes 42/105 misses (40%); also flags 41/215 hits (would be at risk only if a real model overrides).
- **63/105 misses are HIGH confidence (FG confidently wrong) → NOT reachable by the gate**:
  NO-ANSWER 39 (e.g. XO→kinase via pyrimidine), FALSE-POS 24 (e.g. COMT loses to CA/NR at score 6+).
- Implication: a confidence gate alone cannot trigger the FALSE-POSITIVE re-rank cases (the wrong
  winner is *strong*). Reaching the 63 confident-wrong misses needs a parallel 3D path + meta-reconciler
  that may override high-confidence FG calls — which carries regression risk and must be validated hard.
  Method map (from failure-mode analysis): shape/ProLIF = retrieval for NO-ANSWER; Gnina rescore =
  re-rank for FALSE-POS (but needs an always-on or top-K trigger, not the confidence gate).

## ProLIF fallback — runs end-to-end (`utils/fallback_3d.py`, 2026-06-18)

The old OpenBabel protein-prep blocker is RESOLVED and the serine-protease reference library
is BUILT. Pipeline (`ProLIFFallback`, all implemented, not stubbed):
`_embed_3d` (**pH 7.4 protonation via `_protonate_smiles` — same OpenBabel model the builder
uses; amidine/guanidine → cationic, COOH → anionic** — then RDKit ETKDGv3+MMFF) →
`_best_match`, which docks the query into **EACH reference's own receptor** (`_dock_and_ifp`:
smina subprocess, `--autobox_ligand`, `--seed 0`; ProLIF Fingerprint vs that receptor) and
takes the max Jaccard. `build_override()` re-rank/merge contract tested.

- **Per-reference docking, NOT a single shared receptor**: each reference IFP lives in its own
  crystal frame (trypsin/thrombin/FXa). Docking the query into one fixed receptor put
  non-trypsin peptidomimetics in the wrong frame — rivaroxaban self-matched at only **0.33**.
  After per-reference docking (redock into 2W26 FXa) rivaroxaban self-matches at **1.0**;
  benzamidine→3PTB **1.0**; ibuprofen (negative) **0.25** → no false positive. Cost: N
  docks/query (~50s each; ~20 min for the 9-ref SP library). Distinct receptors docked once each.
- **IFP keys are CHAIN-STRIPPED** (`ifp_keys_from_fingerprint` → `rsplit('.',1)[0]`): SP family
  shares chymotrypsin numbering (ASP189/SER195/GLY216/219/HIS57 align across trypsin/thrombin/
  FXa), so dropping the crystal chain id makes references cross-PDB comparable. Key form =
  **"RESNAMERESNUM.Interaction"** (e.g. `ASP189.Cationic`) shared by builder and fallback — NOT
  raw bitvectors (unaligned across mols).
- `propose()` swallows all errors → None; **strict `sim > threshold`** (a boundary
  `sim==threshold` maps to score 0.0, a useless proposal — caffeine docks the trypsin S1 at
  ≈0.60). Heavy deps lazy.
- **Still NOT auto-registered** (default hook = no-op stub) → predict/benchmark byte-identical
  (core 190/220 unchanged). Activate: `register_fallback_3d(ProLIFFallback())`. Regression
  surface = the 41 low-confidence HITS (overrides must not break them).
- **`_find_backend()`**: smina 2020.12.10 lives at `<env>/Library/bin/smina.exe`; `shutil.which`
  misses it without `conda activate`, so `_find_backend` also probes
  `sys.prefix`/{Library/bin,bin,Scripts}.

### ⚠ EMPIRICAL RESULT (2026-06-18) — recovers 0/7 serine-protease misses → NO accuracy gain

Ran all 7 SP misses through per-reference docking vs the 9-co-crystal library: best IFP Jaccard
only 0.33–0.60 (CHEMBL323583 0.60, two at 0.50, …), all below the 0.6 fire threshold; 1
(CHEMBL103874) is high-conf wrong (CA 6.14) so the gate never consults the fallback. Even firing
wouldn't flip them — overturning the standing FG top-1 needs sim≈0.72–0.82. The positive controls
(benzamidine/rivaroxaban →1.0) score high only because they ARE in the reference set
(self-docking); the real misses are different chemotypes that don't reproduce a reference binding
mode. **This confirms structurally that the 7 SP misses are real structural misses, not a
missing-pattern gap.** Threshold-lowering won't help (scores too low + would fire on negatives).
Payoff needs a far broader reference library and/or docking into the TRUE target receptor (not a
family proxy), not tuning. **Do NOT register it expecting a benchmark gain.**

| SP miss | FG top-1 (score) | Confidence | Best SP IFP sim | Recovered? |
|---|---|---|---|---|
| CHEMBL103874 | carbonic anhydrase (6.14) | high | — (gate never consults) | ❌ |
| CHEMBL92615 | GPCR (4.30) | low | 0.50 | ❌ |
| CHEMBL323583 | cytochrome P450 (2.50) | low | 0.60 | ❌ |
| CHEMBL4108739 | nuclear receptor (2.38) | low | 0.33 | ❌ |
| CHEMBL285285 | nuclear receptor (2.38) | low | 0.40 | ❌ |
| CHEMBL1682691 | nuclear receptor (2.38) | low | 0.50 | ❌ |
| CHEMBL291026 | nuclear receptor (2.38) | low | 0.50 | ❌ |

## Reference-library builder (`utils/build_prolif_reference.py check|build`)

Builds `db/prolif_reference_ifp.json` from PDB **co-crystals** (real poses → ProLIF IFP; no
docking on the reference side). **Current library = 9 serine-protease co-crystals**: 3PTB/BEN,
1OYT/FSN, 1DWD/MID, 2ZFF/53U, 1F0R/815 + 4 non-benzamidine peptidomimetic FXa
(**2W26/rivaroxaban, 2P16/GG2, 2RA0/JNJ, 2J34/GS6**) for the SP misses with no Arg-mimetic S1 anchor.

- **PDBFixer protein prep is crash-resilient**: `_pdbfixer_worker` runs in a **subprocess** (so a
  native `addMissingAtoms()` segfault on a pathological structure — e.g. FXa 2W26 — can't kill the
  batch); the parent retries WITHOUT addMissingAtoms (skipped atoms are surface side chains, not
  the S1/triad). This recovered 1OYT and 2W26. Ligand still protonated standalone via OpenBabel
  (`_protonate_ligand`).
- **Env (2026-06-18)**: rdkit 2023.09.6 + ProLIF 2.1.0 + MDAnalysis 2.9.0 + pdbfixer + openmm +
  requests + **smina 2020.12.10** all installed. Interactive shell can't `conda activate` — run via
  the full path `C:\Users\User\miniconda3\envs\chem_target\python.exe`. db/prolif_pdb/ +
  db/prolif_reference_ifp.json gitignored (rebuild offline:
  `python utils/build_prolif_reference.py build --target "serine protease"`).

## Tune list (now that it runs)

- **Per-reference docking cost**: N docks/query is the dominant cost. Consider top-K reference
  receptors, parallelising docks, or a cheaper pre-filter (shape) before docking.
- **Receptor prep**: smina is fed a PDBFixer-prepped PDB (no PDBQT). For accuracy add meeko /
  `prepare_receptor` → PDBQT (charges, rigid/flex). Same for the autobox ligand.
- **Confidence gate thresholds** (`target_predictor`): `_CONF_MIN_TOP1_SCORE=3.0`,
  `_CONF_MIN_MARGIN=0.75` — routes 40% of misses + flags 41 hits; tune vs the hit/miss cross-tab.
- **ProLIF similarity→score** (`ProLIFFallback`): `sim_threshold=0.6`, `score_scale=8.0` —
  FALSE-POS recovery needs sim≈0.9+ to beat score-6 winners; NO-ANSWER ≈0.7. Known limit: a query
  docked into a non-native receptor may not reproduce the crystal pose, so some peptidomimetics
  still land < 0.6 (lowering the threshold risks false positives — see ibuprofen 0.25).
- **Reference set coverage**: 9 SP co-crystals; add more non-benzamidine peptidomimetic complexes
  so the remaining SP misses have a near neighbour.
- **Docking determinism**: smina `--seed 0` set; consider `--num_modes>1` + best-IFP-over-poses,
  tune `--exhaustiveness`/`--autobox_add`.
- **Override policy**: `build_override` boosts/inserts one proposal; for FALSE-POS re-rank at scale,
  consider top-K docking + a meta-reconciler (Option 2) rather than the single low-conf gate.
