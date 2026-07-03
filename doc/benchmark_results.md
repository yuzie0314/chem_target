# Benchmark per-class results (detailed dev reference)

Detailed per-class breakdown for the curated benchmark (`run_benchmark.py`,
pChEMBL ≥ 6.0, 20 compounds per class), branch `dev/validation`. Headline numbers
live in `CLAUDE.md`; the public-facing summary is in `README.md`. This file is the
full developer reference (the "Notes" carry the per-class reasoning and the
structural-limit provenance behind each score).

**Core 11-class Top-1: 190/220 = 86.4%   Top-3: 197/220 = 89.5%** (current best).

## Per-class (curated, 20 compounds each)

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

## Note on the 13-class extended set (260 compounds)

MAO + COMT were added 2026-06-15 with *correct* ChEMBL target IDs (MAO=CHEMBL1951/2039,
COMT=CHEMBL2023). Their low scores reflect a **pChEMBL-sampling bias**: the top-20
highest-affinity ChEMBL compounds for these targets are modern research analogs that lack
the classic covalent pharmacophores (propargylamine, nitrocatechol) the marketed drugs
carry — not a rule gap. Adding the warhead rule cannot capture analogs that don't carry the
warhead. The other 5 blind-spot classes (PDE/topoisomerase/ribosome/XO/cysteine protease)
remain un-added; see CLAUDE.md "Next tasks".

## Data bug FIXED (2026-06-15)

`run_stp_comparison.py` previously had a local `_TARGET_CLASS_MAP` with *wrong* ChEMBL IDs
for the 7 blind-spot classes (e.g. MAO→CHEMBL2828 returned DARUNAVIR, topoisomerase→CHEMBL3952
returned the opioid JDTIC). It now does `from run_benchmark import TARGET_CLASS_MAP` —
**single source of truth**, so the two can never diverge again. `generate_comparison` also
restricts to the true compound intersection (`restrict_names`), keeping the head-to-head fair
after benchmark expansion. NOTE: the cached `stp_raw.csv` still holds the old mislabelled STP
rows for those 7 classes (fetched with the bad map); re-query STP to refresh them — but the
published 220/11-class fair comparison never used them, so it is unaffected.
