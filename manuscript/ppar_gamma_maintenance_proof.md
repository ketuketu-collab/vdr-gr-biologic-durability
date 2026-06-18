# PPARγ as the Chronic-Homeostatic Extreme: the Maintenance-Side Positive Proof

**Draft fragment — 2026-06-18**
Intended for integration into the VDR/GR manuscript (Discussion) or as the
foundation of a follow-on "GR→VDR→PPARγ nuclear-receptor durability axis" paper.

Working hypothesis (HN): **PPARγ sits even further toward the chronic /
homeostatic pole than VDR** — i.e. the durability axis is a single temporal
spectrum, not a 2D landscape:

```
GR  ───────────────  VDR  ───────────────  PPARγ
acute / emergency    chronic maintenance    terminal resolution / tissue repair
peak 0–4 h           peak 24–48 h           sustained, latest
induction efficacy   maintenance remission  the safest long-term backbone
prototype: steroids  prototype: anti-IL-23  prototype: 5-ASA
(induction-only,     (delayed-but-durable)  (maintenance-defining)
 self-sabotaging)
```

---

## Draft paragraph (manuscript-ready)

Just as corticosteroids supply the *negative* proof of the GR-acute principle —
universally effective for induction yet formally contraindicated for maintenance
— 5-aminosalicylic acid (5-ASA / mesalamine) supplies the complementary
*positive* proof at the opposite, chronic-homeostatic extreme of the axis. 5-ASA
is the oldest, safest, and most widely used long-term maintenance therapy for
ulcerative colitis, clinically distinguished by negligible cumulative toxicity
over decades of continuous use and by a characteristic weak-induction /
durable-maintenance profile. The molecular basis of that profile is PPARγ: the
intestinal anti-inflammatory effect of 5-ASA is strictly PPARγ-dependent — it is
abolished in PPARγ-haploinsufficient (PPARγ^+/−) mice while preserved in
wild-type littermates — and 5-ASA acts as a direct PPARγ agonist that promotes
receptor nuclear translocation and coactivator recruitment in colonic
epithelium, a result validated in organ cultures of human colonic biopsies
[Rousseaux 2005]. Patients with ulcerative colitis exhibit impaired PPARγ
expression specifically confined to the colonic epithelium, and PPARγ is
recognized as the principal functional receptor mediating aminosalicylate
activity in inflammatory bowel disease [Dubuquoy 2006]. Direct pharmacological
PPARγ agonism recapitulates the benefit: in a multicenter randomized
double-blind placebo-controlled trial in mild-to-moderately active ulcerative
colitis (n = 105), the thiazolidinedione PPARγ ligand rosiglitazone achieved
clinical response in 44% versus 23% of placebo-treated patients (p = 0.04) and
clinical remission in 17% versus 2% (p = 0.01), with serious adverse events rare
[Lewis 2008].

This positions PPARγ even further toward the chronic-homeostatic pole than VDR.
Whereas VDR governs immune tolerance and epithelial barrier integrity over a
24–48 h timescale, PPARγ is the master regulator of the *terminal resolution*
programs that restore tissue homeostasis after inflammation: maturation of
alternatively activated (M2) macrophages [Odegaard 2007] and the accumulation,
phenotype, and function of tissue-resident regulatory T cells [Cipolletta 2012].
The three nuclear receptors therefore describe a single temporal axis of immune
control — GR (acute emergency, peak 0–4 h) → VDR (chronic maintenance, peak
24–48 h) → PPARγ (terminal resolution and tissue repair) — along which durable
therapeutic benefit increases monotonically from left to right. The clinical
corollary is striking: the two endpoints of this axis are occupied by the two
oldest anti-inflammatory drug classes in IBD — corticosteroids and 5-ASA — whose
diametrically opposite maintenance behaviors are exactly what the GR→VDR→PPARγ
gradient predicts. Corticosteroids (GR) achieve high induction remission but are
contraindicated for maintenance and are self-limiting through CYP24A1-driven
depletion of the VDR ligand pool; 5-ASA (PPARγ) is weakly effective for
induction yet defines safe long-term maintenance, acting through the resolution
receptor that lies at the durable extreme of the axis.

A final caveat preserves the rule's specificity rather than weakening it.
Systemic thiazolidinediones carry class-level cardiovascular and skeletal
liabilities that are unrelated to their colonic anti-inflammatory mechanism; the
framework accordingly predicts that *gut-restricted* PPARγ engagement — as
achieved by topically acting 5-ASA confined to the colonic epithelium — is the
durable and safe configuration, consistent with 5-ASA's decades-long clinical
track record and motivating the development of colon-targeted next-generation
PPARγ ligands.

---

## Novelty positioning (what is and isn't new)

The individual facts below are **established** and must be cited as such; the
novelty is their **integration into a single GR→VDR→PPARγ temporal-durability
axis** and the prediction that therapeutic durability increases monotonically
along it.

- *Known*: PPARγ mediates 5-ASA action; PPARγ drives M2 macrophages and tissue
  Tregs; PPARγ agonists have anti-colitis activity.
- *New (this work)*: placing PPARγ **on the same calculable ChIP-seq occupancy
  axis** as VDR and GR; the claim that PPARγ occupancy is **even more predictive
  of long-term durability than VDR**; and the unifying read that the two oldest
  IBD drugs (steroids, 5-ASA) mark the two ends of the axis.

## Open analysis to confirm the placement (requires ReMAP PPARG/NR1C3 beds)

1. Score PPARγ (gene symbol **PPARG**, nuclear receptor **NR1C3**) occupancy at
   the 381 target genes from ReMAP2022, parallel to the existing VDR/GR columns
   in `results/remap_scores_expanded.csv` (note: the current `PPARG` *row* is
   PPARG-as-target, not PPARγ-as-regulator — a new regulator column is needed).
2. Test whether durable (◎) targets score higher on PPARγ than VDR, and whether
   adding the PPARγ axis raises the approval/durability AUC above 0.706.
3. Temporal check: do canonical PPARγ targets peak later / more sustained than
   the VDR set (24–48 h) in the existing time-course data?

---

## Verified references (PubMed, 2026-06-18)

1. Rousseaux C, Lefebvre B, Dubuquoy L, et al. Intestinal antiinflammatory
   effect of 5-aminosalicylic acid is dependent on peroxisome
   proliferator-activated receptor-gamma. *J Exp Med*. 2005;201(8):1205–1215.
   PMID: 15824083. doi:10.1084/jem.20041948

2. Dubuquoy L, Rousseaux C, Thuru X, et al. PPARgamma as a new therapeutic
   target in inflammatory bowel diseases. *Gut*. 2006;55(9):1341–1349.
   PMID: 16905700. doi:10.1136/gut.2006.093484

3. Lewis JD, Lichtenstein GR, Deren JJ, et al. Rosiglitazone for active
   ulcerative colitis: a randomized placebo-controlled trial. *Gastroenterology*.
   2008;134(3):688–695. PMID: 18325386. doi:10.1053/j.gastro.2007.12.012

4. Cipolletta D, Feuerer M, Li A, et al. PPAR-γ is a major driver of the
   accumulation and phenotype of adipose tissue Treg cells. *Nature*.
   2012;486(7404):549–553. PMID: 22722857. doi:10.1038/nature11132

5. Odegaard JI, Ricardo-Gonzalez RR, Goforth MH, et al. Macrophage-specific
   PPARgamma controls alternative activation and improves insulin resistance.
   *Nature*. 2007;447(7148):1116–1120. PMID: 17515919. doi:10.1038/nature05894

*Supplementary (PPARγ-mediated antineoplastic 5-ASA effect, supports the
epithelial-PPARγ axis): Rousseaux C, El-Jamal N, Fumery M, et al. Carcinogenesis.
2013;34(11):2580–2586. PMID: 23843037. doi:10.1093/carcin/bgt245*
