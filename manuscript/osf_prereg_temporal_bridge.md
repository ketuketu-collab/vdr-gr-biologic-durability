# OSF Preregistration Draft
# Temporal Ordering of Nuclear-Receptor Programs: PPARγ as the Induction→Maintenance Bridge

**Date**: 2026-06-19 (to be timestamped after the Phase-0 pilot locks conditions)
**Investigator**: Hiroyuki Nagashima, Sapporo Medical University
**Companion analysis (frozen)**: `scripts/16_temporal_bridge.py`; protocol
`manuscript/wet_protocol_temporal_bridge.md`

---

## Title
Glucocorticoid, PPARγ and Vitamin-D Receptor Target Programs Fire in Temporal
Order in Human Macrophages: PPARγ Marks the Resolution-Phase Bridge Between
Induction and Maintenance.

---

## Background & rationale (already established; not under test here)
From public ChIP-seq (ReMAP2022) we found PPARγ (NR1C3) occupancy: (i) predicts
long-term durability of biologic targets as well as or better than VDR (durability-
tier separation p = 0.023; IBD UC+CD maintenance-remission Spearman r = 0.899),
(ii) does NOT predict regulatory approval (AUC = 0.38), and (iii) is positioned
equidistant from GR and VDR in occupancy space (PPARγ–GR 0.34, PPARγ–VDR 0.36),
with its durability signal independent of GR (partial Spearman 0.84). These
occupancy results are exploratory/observational and are NOT the object of this
preregistration. This study prospectively tests the **temporal** prediction they
motivate.

---

## Hypotheses

### Primary hypothesis (confirmatory)
In a human-macrophage inflammatory time course, the expression peak time of
ligand-driven target-gene modules follows the order **GR < PPARγ < VDR**:
- GR-module targets peak early (acute; 0–4 h)
- PPARγ-module targets peak in the middle (resolution/bridge; ~8–12 h)
- VDR-module targets peak late (maintenance; 24–48 h)

### Directional sub-hypotheses
- H1: PPARγ-module median peak time is **later** than GR-module.
- H2: PPARγ-module median peak time is **earlier** than VDR-module.

---

## Design (confirmatory run; conditions locked by an exploratory Phase-0 pilot)
THP-1-derived macrophages, single-ligand arms on a constant LPS inflammatory
background (no ligand co-administration; LPS is a fixed backdrop, not a crossed
factor). Arms: (1) LPS+vehicle, (2) LPS+dexamethasone, (3) LPS→IL-4 at 6 h,
(4) LPS+vitamin D. Timepoints 0/2/4/8/12/24/48 h (final set may be refined in
the pilot and will be fixed before sequencing). n = 3 biological replicates.
Read-out: 3'-tag RNA-seq (gene-level counts). Translational confirmation:
RT-qPCR of key genes in primary MDM (≥3 donors).

*Exploratory (NOT preregistered):* the Phase-0 pilot (condition optimization:
DEX dose, 25-D3 vs 1,25-D3, resolution trigger, timepoint placement) and any
combined-cocktail "integration" course.

---

## Variables

### Module gene lists (FROZEN; as encoded in scripts/16_temporal_bridge.py)
- **GR (acute):** TSC22D3, FKBP5, DUSP1, ZBTB16, PER1, KLF13, TXNIP, DDIT4
- **PPARγ (bridge):** CD36, MRC1, ANGPTL4, FABP4, MERTK, ALOX15, CD163, PPARG, LPL, PLIN2
- **VDR (maintenance):** CYP24A1, CAMP, NOD2, TLR10, DEFB4A, IL37, CD14

### Outcome
Per-gene expression peak time = hour of maximum expression (argmax over
timepoints), in the gene's own ligand arm; aggregated per module (median).

---

## Analysis plan (FROZEN)
1. 3'-tag reads → UMI-collapsed gene counts → log2FC vs 0 h.
2. Run `scripts/16_temporal_bridge.py` on the gene × timepoint matrix.
3. Pre-specified tests:
   - Directional Mann-Whitney: PPARγ peak times > GR (H1); PPARγ < VDR (H2).
   - Kruskal-Wallis across the three module peak-time distributions.
4. (Secondary) continuous peak-time + CI via impulse/spline fit.

---

## Inference criteria (success)
The primary hypothesis is **supported** iff ALL hold:
- median peak order GR < PPARγ < VDR, AND
- H1 (PPARγ later than GR) p < 0.05, AND
- H2 (PPARγ earlier than VDR) p < 0.05, AND
- Kruskal-Wallis across modules p < 0.05.

Partial support (e.g. PPARγ < VDR only) will be reported as such, not as
confirmation.

---

## Sample size / power
4 arms × 7 timepoints × 3 replicates = 84 RNA-seq samples. Module-level peak-time
ordering aggregates 7–12 genes per module; this design is powered for the
ordering test (pilot confirms effect direction before commitment). No interim
analyses; no optional stopping.

---

## What would falsify the hypothesis
- PPARγ-module peak ≤ GR-module peak (no bridge; PPARγ co-temporal with acute), or
- PPARγ-module peak ≥ VDR-module peak (PPARγ not earlier than maintenance), or
- no significant peak-time separation across modules.

---

## Caveats (pre-stated)
- THP-1 M2/resolution (PPARγ) may be blunted vs primary; MDM qPCR confirms.
- Module gene lists are frozen here to prevent post-hoc selection.
- Cross-stimulus peak comparison is mitigated by the shared LPS background and
  the within-arm ±ligand contrast.

---

## Timeline
Phase-0 pilot (exploratory) → freeze conditions → **timestamp this preregistration**
→ main RNA-seq → analysis with the frozen `scripts/16`.
