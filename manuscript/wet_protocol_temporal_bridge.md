# Wet Protocol — Temporal Bridge of Nuclear-Receptor Programs (GR → PPARγ → VDR)

**Draft design v1 — 2026-06-19 (HN)**
Reusable experimental template. A single stimulus×time-course design that, once
fixed, transfers to organoid / patient-biopsy / mouse-KO arms (each adds one
figure at marginal cost). Companion analysis: `scripts/16_temporal_bridge.py`.

---

## 1. Hypothesis & success criterion

**Hypothesis.** During an inflammation→resolution course, nuclear-receptor target
programs fire in temporal order: GR (acute) → PPARγ (resolution/bridge) → VDR
(maintenance). PPARγ marks the induction→maintenance hand-off.

**Primary readout.** Per-gene expression peak time, aggregated per module.
**Success.** median peak hour GR < PPARγ < VDR, with PPARγ significantly later
than GR and earlier than VDR (directional Mann-Whitney; Kruskal-Wallis across
modules). Target window: GR 0–4 h, PPARγ ~8–12 h, VDR 24–48 h.

---

## 2. The design hinge — making all three programs fire in ONE course

Each receptor needs its ligand, so a naïve LPS course will NOT light up GR/VDR.
The key design choice is to let each ligand be **available or generated in its
natural temporal window**:

**Primary design (physiological resolution course).**
Human macrophages stimulated to inflame then resolve, in medium that permits
endogenous/local ligand generation:
- **Acute/GR:** LPS (or zymosan) ± low-dose dexamethasone (10 nM) at t0 → GR
  immediate-early targets (GILZ/FKBP5/DUSP1) peak 0–4 h.
- **Resolution/PPARγ:** drive resolution by efferocytosis (add apoptotic
  neutrophils/PMN at ~6 h) or IL-4 (20 ng/mL) → endogenous 15-LOX/15d-PGJ₂
  generation activates PPARγ; M2/efferocytosis targets (CD36/MRC1/MERTK/ALOX15)
  peak ~8–12 h.
- **Maintenance/VDR:** supplement **25-hydroxyvitamin D₃ (100 nM)** in the medium;
  activated macrophages express CYP27B1 and convert it locally to calcitriol,
  so VDR targets (CYP24A1/CAMP/NOD2) rise late, 24–48 h. (This local-conversion
  step is genuine macrophage biology and is the elegant part of the design.)

**Validation arms (simpler, stimulus-specific — confirm each module separately).**
- GR module: dexamethasone ± LPS time course (cf. GSE93735).
- VDR module: 1,25-D₃ (or 25-D₃ + activation) time course (cf. GSE189984).
- These establish per-module kinetics under a clean single ligand; the primary
  course establishes the ordering within one biological context.

> Pitfall to state explicitly: if modules are read across different stimuli/cell
> types, ordering is suggestive, not within-stimulus proof. The primary design
> exists precisely to avoid that. Refine exact concentrations/timing with the wet
> collaborator.

---

## 3. Sampling plan

| Axis | Choice |
|---|---|
| Cell model | THP-1-derived macrophages (PMA-differentiated) for MVE → primary human MDM (buffy coat, n≥3 donors) for the paper figure |
| Timepoints | 0, 2, 4, 8, 12, 24, 48 h (9-point 0/1/2/4/6/8/12/24/48 h if peak-time CIs wanted) |
| Replicates | n = 3 (biological); n = 4 for the robust tier |
| Arms | ± resolution stimulus (2) for the primary course |

---

## 4. Sample count & cost (ballpark, JPY; ¥150/$)

Per-sample all-in: 3'-tag/BRB-seq ≈ ¥6,000; full mRNA-seq ≈ ¥20,000.
Reagents (cells, culture, RNA extraction) ≈ ¥150,000.

| Scenario | n | 3'-tag total | full mRNA-seq |
|---|---|---|---|
| A — MVE (1 arm ×7 tp ×3) | 21 | **¥276k ($1.8k)** | ¥570k ($3.8k) |
| B — standard (2 arm ×7 tp ×3) | 42 | **¥402k ($2.7k)** | ¥990k ($6.6k) |
| C — robust (2 arm ×7 tp ×4) | 56 | ¥486k ($3.2k) | ¥1.27M ($8.5k) |

**Platform: 3'-tag RNA-seq (BRB-seq/QuantSeq).** Time courses only need gene-level
counts → 5M reads/sample suffices, cheap multiplexed library prep. Half the cost
of full mRNA-seq with negligible information loss for this question.

**Recommended route:** ① RT-qPCR pilot (~25 genes, ¥5–10万, 2 wk) to confirm the
ordering → ② Scenario B 3'-tag RNA-seq (¥40万) for the publication figure.
**~¥45–50万 (~$3k) buys one temporal-bridge figure.**

---

## 5. Gene modules / qPCR panel (also hard-coded in scripts/16)

- **GR (acute):** TSC22D3 (GILZ), FKBP5, DUSP1, ZBTB16, PER1, KLF13, TXNIP, DDIT4
- **PPARγ (bridge):** CD36, MRC1/CD206, ANGPTL4, FABP4, MERTK, ALOX15, CD163,
  PPARG, LPL, PLIN2 — note ALOX15 doubles as the enzyme making the endogenous
  PPARγ ligand (15-HETE/15d-PGJ₂), a feed-forward marker of the bridge
- **VDR (maintenance):** CYP24A1, CAMP, NOD2, TLR10, DEFB4A, IL37, CD14
- Housekeeping for qPCR: GAPDH, ACTB, B2M (geometric mean)

---

## 6. Analysis & statistics

- Quantify (3'-tag) → gene × timepoint count matrix → CPM/log2FC vs 0 h.
- `python scripts/16_temporal_bridge.py --matrix counts_timecourse.csv` →
  per-gene peak hour, module medians, directional Mann-Whitney, Kruskal-Wallis,
  peak-time figure.
- Optional: model each gene with a smooth/impulse fit (e.g. ImpulseDE2) for a
  continuous peak-time estimate + CI rather than discrete argmax.
- Power: n = 3 × 7 tp is ample for module-level ordering; the peak-time argmax is
  robust to per-gene noise because it aggregates 7–12 genes per module.

---

## 7. Reusability roadmap (the multi-paper value)

The same stimulus×time template, swapping only the biological system, yields a
figure per arm at marginal cost:

1. **Cell line (THP-1)** — MVE / assay development.
2. **Primary human MDM** — main in-vitro figure.
3. **Patient-derived colonic organoids** (IBD vs control) — translational; does
   the bridge operate in epithelium and is it impaired in IBD?
4. **Patient biopsies** (active vs resolving vs stable remission) — does PPARγ
   peak at the resolution/transition phase in vivo? (HN has endoscopy access.)
5. **Mouse macrophage-PPARγ KO** (LysM-Cre × Pparg^fl) — causal: does removing
   the bridge abolish resolution and durable remission?

Arms 1–4 are the dry+wet package for the standalone mechanism paper; arm 5 is the
causal capstone / follow-on.

---

## 8. Key materials (indicative)

LPS (E. coli O111:B4), dexamethasone, recombinant IL-4, 25-OH-vitamin D₃,
PMA (for THP-1), Ficoll/CD14 beads (MDM isolation), apoptotic PMN (efferocytosis),
RNA extraction kit, 3'-tag library kit (Lexogen QuantSeq / BRB-seq), qPCR
mastermix + primers for the panel above.

---

## 9. Pre-registration

Consistent with the lab's OSF practice (osf.io/tnp63), pre-register the primary
hypothesis (GR < PPARγ < VDR median peak ordering), the module gene lists, the
success criterion, and the analysis script before sequencing.
