# Protocol — Temporal Bridge of Nuclear-Receptor Programs (GR → PPARγ → VDR)

**Bench protocol v2 (finalized design) — 2026-06-19, H. Nagashima**
Reusable stimulus×time-course template. Companion analysis: `scripts/16_temporal_bridge.py`.
Hand to wet-lab collaborator / core facility; concentrations are validated
starting points to optimize in pilot.

---

## 1. Objective, hypothesis, endpoint

**Objective.** Test, in one biological context, whether nuclear-receptor target
programs fire in temporal order during an inflammation→resolution→maintenance
course of human macrophages.

**Hypothesis.** GR targets peak early (acute), PPARγ targets in the middle
(resolution/bridge), VDR targets late (maintenance).

**Primary endpoint.** Per-gene expression peak time, aggregated per module.

**Success criterion.** median peak hour **GR < PPARγ < VDR**, PPARγ significantly
later than GR and earlier than VDR (directional Mann-Whitney) + Kruskal-Wallis
across modules p < 0.05. Target windows: GR 0–4 h, PPARγ 8–12 h, VDR 24–48 h.

**Design logic.** All three ligands are supplied at t0 (GR and VDR ligands are
absent in vitro and must be added; the PPARγ ligand is generated endogenously by
the resolution program). Because ligand availability is therefore simultaneous,
an observed GR<PPARγ<VDR ordering reflects the **intrinsic kinetics** of each
program, not staggered ligand delivery — the strongest form of the argument.

---

## 2. Experimental scheme

```
 Macrophages (THP-1 → MDM)              harvest at 0 / 2 / 4 / 8 / 12 / 24 / 48 h
 │
 t0 ── LPS 100 ng/mL                (inflammation)
     + Dexamethasone 100 nM         (GR ligand → GILZ/FKBP5 early)
     + 25-OH-vitamin D3 100 nM      (→ CYP27B1 → calcitriol → CYP24A1 late)
 t6 ── IL-4 20 ng/mL  OR  apoptotic PMN (efferocytosis)
                                    (resolution → endogenous 15-HETE → PPARγ mid)
```

---

## 3. Conditions & sample count

**Main RNA-seq (3'-tag) — combined course, sample-efficient (read all 3 modules):**

| Arm | t0 | t6 | n (×7 tp ×3 rep) |
|---|---|---|---|
| 1 Bridge course | LPS + Dex + 25-D3 | + IL-4 (resolution) | 21 |
| 2 Inflammation control | LPS only | — | 21 |
| **Total** | | | **42** |

**Deconvolution (RT-qPCR only — cheap, confirms module attribution):**

| Single-ligand arm | stimulus | read-out module |
|---|---|---|
| GR | Dexamethasone 100 nM | GILZ/FKBP5/DUSP1 |
| VDR | 1,25-D3 10 nM (or 25-D3 + LPS) | CYP24A1/CAMP/NOD2 |
| PPARγ | Rosiglitazone 1 µM (positive control) | CD36/MRC1/ALOX15 |

Each single-ligand arm: 0/2/4/8/12/24/48 h × n3, qPCR panel only.

---

## 4. Materials

**Cells.** THP-1 (ATCC TIB-202) for MVE; human monocyte-derived macrophages
(MDM) from buffy coat / leukapheresis (Ficoll → CD14+ selection) for the paper.

**Stimuli / ligands.**
- LPS, E. coli O111:B4 (100 ng/mL)
- Dexamethasone (100 nM; pilot 10–100 nM)
- 25-hydroxyvitamin D3 (100 nM) — local CYP27B1 activation; or 1,25-D3 (10 nM)
- Recombinant human IL-4 (20 ng/mL)
- Rosiglitazone (1 µM; PPARγ positive control)
- (efferocytosis option) apoptotic human neutrophils, 5:1 PMN:macrophage
- PMA (25–50 ng/mL, THP-1 differentiation), M-CSF (50 ng/mL, MDM)

**Kits.** RNA extraction column kit; 3' RNA-seq library kit with UMIs (Lexogen
QuantSeq 3' FWD or BRB-seq); qPCR mastermix + primers (panel §6).

---

## 5. Methods

**5.1 Cell preparation.**
- *THP-1:* differentiate with PMA 25–50 ng/mL ×48–72 h; wash; rest 24 h in
  PMA-free medium before stimulation.
- *MDM:* PBMC by Ficoll; CD14+ monocytes (beads); differentiate with M-CSF
  50 ng/mL ×6–7 d. Use ≥3 independent donors as biological replicates.

**5.2 Seeding.** 0.5–1×10⁶ macrophages/well (12-well), rest overnight, serum
conditions held constant across all wells/timepoints.

**5.3 Stimulation.** Apply t0 cocktail per arm (§3). At t6 add resolution trigger
(IL-4 or apoptotic PMN) to Arm 1. Stagger start times so all timepoints harvest
together (recommended) — i.e. plate a separate well per timepoint.

**5.4 Harvest.** At each timepoint aspirate, lyse directly in column-kit lysis
buffer (or TRIzol), snap-freeze, −80 °C. The 0 h sample = pre-stimulation.

**5.5 RNA QC.** Quantify (Qubit), integrity (Bioanalyzer/TapeStation); require
RIN ≥ 8 for 3'-tag input.

**5.6 Library & sequencing.** 3'-tag library with UMIs; pool 42 samples; single-end
75–100 bp; **~5 M reads/sample** (≥3 M usable). NextSeq/NovaSeq.

**5.7 qPCR deconvolution.** Reverse-transcribe; SYBR/TaqMan for panel genes;
normalize to geometric mean of GAPDH/ACTB/B2M; ΔΔCt vs 0 h.

---

## 6. Gene modules / qPCR panel (also in scripts/16)

- **GR (acute):** TSC22D3 (GILZ), FKBP5, DUSP1, ZBTB16, PER1, KLF13, TXNIP, DDIT4
- **PPARγ (bridge):** CD36, MRC1/CD206, ANGPTL4, FABP4, MERTK, ALOX15, CD163,
  PPARG, LPL, PLIN2  *(ALOX15 = enzyme making the endogenous PPARγ ligand → bridge feed-forward marker)*
- **VDR (maintenance):** CYP24A1, CAMP, NOD2, TLR10, DEFB4A, IL37, CD14
- **Reference:** GAPDH, ACTB, B2M

---

## 7. Analysis & statistics

1. 3'-tag reads → UMI-collapsed gene counts → CPM, log2FC vs 0 h.
2. `python scripts/16_temporal_bridge.py --matrix counts_timecourse.csv` →
   per-gene peak hour, module medians, directional Mann-Whitney, Kruskal-Wallis,
   peak-time figure.
3. Optional continuous peak-time + CI via impulse/spline fit (e.g. ImpulseDE2).
4. **Power.** 7 timepoints × n3 is ample for module-level ordering; argmax
   peak-time aggregates 7–12 genes/module, robust to per-gene noise.

---

## 8. Timeline (indicative)

| Phase | Weeks |
|---|---|
| qPCR pilot (THP-1, ordering sanity) | 1–3 |
| Main course + harvest (MDM, ≥3 donors) | 4–7 |
| RNA QC + library + sequencing | 8–10 |
| Analysis + figure | 11–12 |

---

## 9. Budget (ballpark, JPY; ¥150/$)

Per-sample all-in: 3'-tag ≈ ¥6,000; full mRNA-seq ≈ ¥20,000. Reagents ≈ ¥150k.

| Scenario | n | 3'-tag | full |
|---|---|---|---|
| MVE (1 arm) | 21 | ¥276k | ¥570k |
| **Standard (this protocol, 2 arm)** | **42** | **¥402k (~$2.7k)** | ¥990k |
| Robust (2 arm ×4 rep) | 56 | ¥486k | ¥1.27M |

Platform: **3'-tag/BRB-seq** — gene-level counts suffice for peak timing, ½ the
cost of full mRNA-seq, all panel genes (incl. CYP24A1) fully quantified.
qPCR deconvolution arms ≈ ¥5–10万.

---

## 10. Risks & optimization (pilot first)

- **Dex blunts inflammation/resolution.** Dex is anti-inflammatory; high doses may
  suppress LPS/IL-4 programs. Pilot 10 vs 100 nM; consider transient Dex (wash at
  4 h) or a Dex-as-separate-arm design if the combined course is dominated.
- **VDR ligand.** If macrophage CYP27B1 conversion of 25-D3 is weak in vitro, use
  1,25-D3 (10 nM) directly for the VDR arm.
- **Resolution trigger.** IL-4 is simpler/reproducible; efferocytosis is more
  physiological but variable — start with IL-4.
- **3' bias / RNA quality.** Require RIN ≥ 8; for future FFPE biopsy arms switch to
  an FFPE-compatible 3' kit or RNA-capture.
- **Donor variability (MDM).** Use ≥3 donors; block by donor in the model.

---

## 11. Reusability roadmap (multi-paper value)

Same stimulus×time template, swap the biological system → one figure/arm at
marginal cost:
1. THP-1 — MVE/assay dev.
2. Primary MDM — main in-vitro figure.
3. IBD vs control colonic **organoids** — epithelial bridge; impaired in IBD?
4. Patient **biopsies** (active / resolving / stable remission) — PPARγ peaks at
   the transition in vivo? (endoscopy access.)
5. Macrophage-**PPARγ KO mouse** (LysM-Cre × Pparg^fl) — causal: remove the bridge
   → resolution + durable remission fail.
Arms 1–4 = dry+wet package for the standalone mechanism paper; arm 5 = causal capstone.

---

## 12. Pre-registration

Per lab OSF practice (osf.io/tnp63): pre-register the primary hypothesis
(GR < PPARγ < VDR median peak ordering), module gene lists, success criterion,
and `scripts/16` analysis BEFORE sequencing.
