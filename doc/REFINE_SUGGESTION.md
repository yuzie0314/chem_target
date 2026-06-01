# Reverse Virtual Screening (RVS) Bottleneck Analysis

## Current Architecture

Current RVS framework:

* Pharmacophore matching
* Interaction DB matching
* 3D binding state inference

Core issue:

```text
Ligand → Target mapping is many-to-many
```

Certain target families inherently share highly similar ligand interaction patterns.

---

# 1. Structural Ceiling Problems (Hard to Solve)

These are not mainly engineering failures.

They are caused by:

```text
Low information entropy in ligand space
```

---

# mTOR / Kinase Family

## Problem

Most ligands are ATP-competitive inhibitors.

Shared features across kinases:

* Hinge binder
* HBA/HBD pair
* Flat aromatic scaffold
* Hydrophobic cap

After abstraction:

```text
[Aromatic]
[HBA]
[Hydrophobe]
```

becomes nearly identical across:

* mTOR
* PI3K
* EGFR
* CDK
* VEGFR

---

## Why Difficult

Ligands are intentionally designed to overlap across kinases.

Target prediction becomes:

```text
Non-identifiable
```

---

## Recommended Direction

### 1. Add Interaction Geometry

Do not only use:

```text
Feature presence
```

Add:

* Hinge angle
* Donor/acceptor vector
* Aromatic plane distance
* Solvent exposure
* Spatial tolerance

Geometry carries more information than FG presence.

---

### 2. Hierarchical Target Scoring

Avoid forced single-target classification.

Use:

```text
Kinase family
→ subfamily reranking
```

Example:

```text
Kinase: 0.82
 ├─ PI3K/mTOR: 0.31
 ├─ EGFR: 0.28
 ├─ CDK: 0.21
```

---

### 3. Add Shape/Flexibility Features

Useful descriptors:

* PMI
* Shape anisotropy
* Rotatable bonds
* Conformational entropy

Kinase inhibitors are often:

* Flat
* Rigid
* Hinge-centric

GPCR ligands are usually:

* Flexible
* Globular

---

# Adenosine Receptor

## Problem

Many ligands are non-purine scaffolds.

Shared interaction topology:

* π-π stacking
* Hydrophobic cavity
* One polar anchor

FG-based discrimination becomes weak.

---

## Recommended Direction

Do not focus on SMARTS expansion.

Instead add:

* Aromatic centroid network
* Electrostatic surface moments
* Hydrophobe spacing
* Charge distribution

Focus on:

```text
Interaction topology
```

rather than FG identity.

---

# 2. Engineering / Scoring Problems (More Solvable)

These are mainly scoring bias issues.

---

# HDAC vs GPCR Competition

## Problem

Motif:

```text
Tertiary amine + phenyl ring
```

strongly activates GPCR scoring.

HDAC signal gets suppressed.

---

## Root Cause

Current scoring likely resembles:

```text
Σ(IDF-weighted pharmacophore votes)
```

Problem:

Biologically important motifs may have low IDF.

Example:

```text
Hydroxamate
```

is common in HDAC ligands,
therefore IDF becomes weak,
despite being mechanistically critical.

---

## Recommended Direction

### 1. Separate Importance From Rarity

Do NOT assume:

```text
importance = rarity
```

Use:

```text
final_score =
importance_weight × rarity_weight
```

---

### 2. Mechanistic Interaction Bonus

Boost high-information interactions:

* Zn chelation
* Bidentate H-bond
* Catalytic interaction
* Salt bridge pair

These should outweigh generic aromatic overlap.

---

### 3. Saturating GPCR Motif Scoring

Prevent:

```text
amine + aromatic
```

from dominating all predictions.

Use diminishing returns:

```text
First motif → strong gain
Second motif → smaller gain
```

---

# Serine Protease vs GPCR

## Problem

Shared motifs:

* Amine
* Carboxylate

cause GPCR dominance.

However true protease recognition depends more on:

```text
Peptide-like geometry
```

than FG presence.

---

## Recommended Direction

### 1. Pseudo-Backbone Detection

Detect:

* Amide repetition
* Directionality
* Distance matrix patterns
* Backbone vectors

rather than isolated SMARTS.

---

### 2. Shape Grammar

Differentiate:

```text
Linear peptide-like
vs
Compact GPCR ligand
```

Useful descriptors:

* Graph diameter
* Radius of gyration
* PMI
* Asphericity

---

### 3. Context-Aware SMARTS

Do not detect:

```text
Amine
```

alone.

Instead detect:

* Amide-adjacent amine
* Beta-carbonyl amine
* Peptidomimetic context

Context matters more than isolated FG.

---

# Highest ROI Improvements

## Tier 1 (Most Important)

### Interaction Geometry

Most impactful upgrade.

---

### Feature Context

Move beyond FG presence.

---

### Hierarchical Target Scoring

Target family → subtype ranking.

---

### Mechanistic Interaction Weighting

Examples:

* Zn chelation
* Hinge binding
* Catalytic triad interaction

---

# Tier 2

### Shape Descriptors

* PMI
* Globularity
* Asphericity

---

### Ligand Flexibility

* Rotatable bonds
* Entropy

---

# Tier 3

### More SMARTS Expansion

Lowest ROI.

Main limitation is:

```text
Contextual interaction representation
```

not lack of FG dictionaries.

---

# Suggested System Redefinition

Current system:

```text
Ligand FG matcher
```

Target architecture:

```text
Interaction field inference engine
```

Core components:

* FG
* Geometry
* Topology
* Shape
* Mechanistic priors
* Hierarchical target modeling
* Interaction context
