# Protocol — Temporal Bridge of Nuclear-Receptor Programs (GR → PPARγ → VDR)

**Bench protocol v3 (finalized design) — 2026-06-19, H. Nagashima**
Reusable stimulus×time-course template. Companion analysis: `scripts/16_temporal_bridge.py`.
Hand to wet-lab collaborator / core facility; concentrations are validated
starting points, locked by the Phase-0 pilot (§2).

---

## 1. Objective, hypothesis, endpoint

**Objective.** Test whether nuclear-receptor target programs fire in temporal
order on a common inflammatory background in human macrophages.

**Hypothesis.** GR targets peak early (acute), PPARγ targets in the middle
(resolution/bridge), VDR targets late (maintenance).

**Primary endpoint.** Per-gene expression peak time, aggregated per module.

**Success criterion.** median peak hour **GR < PPARγ < VDR**, PPARγ significantly
later than GR and earlier than VDR (directional Mann-Whitney) + Kruskal-Wallis
across modules p < 0.05. Target windows: GR 0–4 h, PPARγ 8–12 h, VDR 24–48 h.

**Design logic (corrected).** Each receptor is engaged by its OWN ligand in a
SEPARATE arm (no co-administration — a DEX+VitD+IL-4 cocktail confounds attribution
and DEX would suppress the resolution program). All arms share a constant LPS
inflammatory background (PPARγ resolution requires prior inflammation; IBD is an
inflamed state), so LPS is a fixed backdrop, NOT a crossed factor. Each ligand arm
is read against the LPS-only arm (the ±ligand contrast). Ordering then reflects the
intrinsic kinetics of each program. LPS-independence of GR/VDR kinetics is
established in the literature and is cited rather than re-proven (optional no-LPS
qPCR spot-check only if a reviewer insists).

---

## 2. Phase 0 — pilot (do this FIRST; ~¥5–10万, 2–3 wk)

De-risk the ¥66万 RNA-seq run. THP-1 macrophages (cheap, no donor variability),
RT-qPCR only, reduced timepoints **0 / 4 / 8 / 24 / 48 h**, n = 3, minimal panel
(GR: GILZ, FKBP5 · PPARγ: CD36, ALOX15 · VDR: CYP24A1, CAMP · ref: GAPDH, ACTB).

Same 4 arms as the main run (§3). Go/no-go + optimization gates:

| Question | Read-out | Decision |
|---|---|---|
| Does the ordering appear? | GILZ early, CD36/ALOX15 mid, CYP24A1 late | GO to RNA-seq if yes |
| Does DEX blunt resolution? | IL-4 arm CD36/ALOX15 induction ± DEX | pick DEX dose (10 vs 100 nM) / transient DEX |
| Does 25-OH-D3 work (CYP27B1)? | CYP24A1 induction with 25-D3 vs 1,25-D3 | switch to 1,25-D3 (10 nM) if weak |
| Best resolution trigger? | CD36/MERTK/ALOX15 kinetics, IL-4 vs efferocytosis | start IL-4; efferocytosis if needed |
| Timepoints well placed? | are peaks captured / saturating? | refine to 7–9 tp for main run |
| RNA quality / primers OK? | RIN, qPCR efficiencies | lock assay |

Only after the pilot locks conditions do you commit to the main RNA-seq.

---

## 3. Main experiment — conditions & sample count

**Single-ligand arms on a constant LPS background (3'-tag RNA-seq):**

| Arm | t0 | t6 | reads out | n (×7 tp ×3) |
|---|---|---|---|---|
| 1 LPS + vehicle | LPS 100 ng/mL | — | baseline (= ±ligand control) | 21 |
| 2 LPS + DEX | LPS + Dex | — | GR module (GILZ/FKBP5) | 21 |
| 3 LPS → IL-4 | LPS | + IL-4 | PPARγ module (CD36/ALOX15) | 21 |
| 4 LPS + VitD | LPS + 25-D3 (or 1,25-D3) | — | VDR module (CYP24A1/CAMP) | 21 |
| **Total** | | | | **84** |

Timepoints 0 / 2 / 4 / 8 / 12 / 24 / 48 h (refine in pilot). Each ligand arm vs
Arm 1 = the ±ligand contrast on a fixed inflamed background.

*Not included (by design):* a ±LPS factorial (would double to 168 — unnecessary,
LPS is a fixed backdrop) and no-LPS single-ligand RNA-seq arms (GR/VDR ligand-driven
kinetics are established; cite, optional qPCR spot-check).

*Optional integration figure (secondary, confounded — label as such):* a single
LPS+DEX+VitD+IL-4 combined course to show the bridge "in one context"; not used for
the primary ordering claim.

---

## 4. Materials

**Cells.** THP-1 (ATCC TIB-202) — pilot + main RNA-seq; human monocyte-derived
macrophages (MDM) from buffy coat (Ficoll → CD14+ selection) — qPCR confirmation only.

**Stimuli / ligands.**
- LPS, E. coli O111:B4 (100 ng/mL)
- Dexamethasone (100 nM; pilot 10–100 nM)
- 25-hydroxyvitamin D3 (100 nM) — local CYP27B1 activation; or 1,25-D3 (10 nM)
- Recombinant human IL-4 (20 ng/mL)
- Rosiglitazone (1 µM; PPARγ positive control)
- (efferocytosis option) apoptotic human neutrophils, 5:1 PMN:macrophage
- PMA (25–50 ng/mL, THP-1 differentiation), M-CSF (50 ng/mL, MDM)

**RNA extraction.** Column or magnetic-bead total-RNA kit with **on-column/in-line
DNase** (genomic-DNA removal is essential for 3'-tag + qPCR). Macrophages are
RNase-rich → a TRIzol/guanidinium-compatible kit is robust:
- Zymo **Direct-zol** (TRIzol lysate → column; good for macrophages), or
  Qiagen **RNeasy** (Mini, or **96-well plate** for the 84–112-sample run), or
  Macherey-Nagel NucleoSpin RNA.
- High throughput: 96-well plate or magnetic-bead (e.g. Zymo Quick-RNA 96 /
  KingFisher) to process all timepoints in one batch and minimize handling drift.
- Target ~¥400–700/sample → 84 samples ≈ ¥35k–60k (within the ¥150k reagent line);
  add pilot (~40 qPCR samples).
- Require **RIN ≥ 8** (TapeStation/Bioanalyzer) for 3'-tag input.
- *Cost-saver:* some 3'-tag kits (BRB-seq) accept **crude cell lysate**, skipping
  column purification — confirm with the core facility; can remove the extraction
  step for the RNA-seq arm (keep purification for qPCR/QC).

**Other kits.** 3' RNA-seq library kit with UMIs (Lexogen QuantSeq 3' FWD or
BRB-seq); reverse-transcription + qPCR mastermix + primers (panel §6).

---

## 5. Methods

**Cell-model strategy.** **THP-1 is the primary workhorse** for both the pilot and
the main RNA-seq run — cheap, reproducible, no donor variability, and consistent
with the lab's prior THP-1 dataset (GSE135130). All three programs are functional
in THP-1 (GR→GILZ, VDR→CYP24A1 via endogenous CYP27B1, PPARγ→CD36). Primary MDM is
demoted to a focused **confirmation** (qPCR only, key genes, ≥3 donors) to secure
translational credibility without a second full RNA-seq. Caveat: THP-1 M2/resolution
(PPARγ) can be blunted vs primary — the Phase-0 pilot checks this and, if weak,
raises IL-4 or switches to efferocytosis.

**5.1 Cell preparation.**
- *THP-1 (main):* differentiate with PMA 25–50 ng/mL ×48–72 h; wash; **rest 24–48 h**
  in PMA-free medium before stimulation (let PMA/PKC effects subside — critical).
- *MDM (confirmation only):* PBMC by Ficoll; CD14+ monocytes (beads); differentiate
  with M-CSF 50 ng/mL ×6–7 d; ≥3 donors; qPCR of GILZ/CD36/ALOX15/CYP24A1.

**5.2 Seeding.** 0.5–1×10⁶ macrophages/well (12-well), rest overnight, serum lot
held constant across all wells/timepoints/arms.

**5.3 Stimulation.** Apply t0 stimulus per arm (§3): Arm 1 LPS only; Arm 2 LPS+Dex;
Arm 3 LPS then IL-4 (or apoptotic PMN) at t6; Arm 4 LPS+VitD. Plate a separate well
per timepoint; stagger start times so all timepoints harvest together. Match
vehicle (ethanol/DMSO) concentration across all arms (Dex/VitD/rosi stocks).

**5.4 Harvest.** At each timepoint aspirate, lyse directly in column-kit lysis
buffer (or TRIzol), snap-freeze, −80 °C. The 0 h sample = pre-stimulation.

**5.5 RNA QC.** Quantify (Qubit), integrity (Bioanalyzer/TapeStation); require
RIN ≥ 8 for 3'-tag input.

**5.6 Library & sequencing.** 3'-tag library with UMIs; pool 84 samples; single-end
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
| **Phase 0 — qPCR pilot (THP-1, §2): ordering + condition lock** | 1–3 |
| Main course + harvest (THP-1; 4 arms ×7 tp ×3) | 4–6 |
| RNA QC + library + sequencing (84 samples) | 7–9 |
| Analysis + figure | 10–11 |
| MDM confirmation (qPCR, ≥3 donors, key genes) | parallel 6–10 |

---

## 9. Budget (ballpark, JPY; ¥150/$)

Per-sample all-in: 3'-tag ≈ ¥6,000; full mRNA-seq ≈ ¥20,000. Reagents ≈ ¥150k.

| Phase | n | 3'-tag | full |
|---|---|---|---|
| Phase 0 pilot (qPCR only, THP-1) | — | **¥5–10万** | — |
| **Main run — 4 arms (LPS bg) ×7 tp ×3** | **84** | **¥654k (~$4.4k)** | ¥1.83M |
| Robust (×4 rep) | 112 | ¥822k | ¥2.39M |

Platform: **3'-tag/BRB-seq** — gene-level counts suffice for peak timing, ½ the
cost of full mRNA-seq, all panel genes (incl. CYP24A1) fully quantified. A ±LPS
factorial (168) is deliberately avoided. **Total to a publishable temporal figure:
pilot ¥5–10万 → main ¥66万 ≈ ~¥75万 ($5k).**

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
1. **THP-1 — main RNA-seq** (pilot + full ordering; this protocol).
2. **Primary MDM — qPCR confirmation** (key genes, ≥3 donors; translational check).
3. IBD vs control colonic **organoids** — epithelial bridge; impaired in IBD?
4. Patient **biopsies** (active / resolving / stable remission) — PPARγ peaks at
   the transition in vivo? (endoscopy access.)
5. Macrophage-**PPARγ KO mouse** (LysM-Cre × Pparg^fl) — causal: remove the bridge
   → resolution + durable remission fail.
Arms 1–2 = core dry+wet package for the standalone mechanism paper; 3–4 extend
translation; arm 5 = causal capstone.

---

## 12. Pre-registration

Per lab OSF practice (osf.io/tnp63): pre-register the primary hypothesis
(GR < PPARγ < VDR median peak ordering), module gene lists, success criterion,
and `scripts/16` analysis BEFORE sequencing.
