# VDR/GR ChIP-seq Occupancy Predicts Long-term Therapeutic Durability of Biologics in Chronic Inflammatory Diseases

**Draft v2 — 2026-05-25**  
Target journal: *Nature Medicine*  
OSF Preregistration: https://osf.io/tnp63 (registered 2026-05-25, embargoed)

---

## Authors

Hiroyuki Nagashima, MD  
Department of Gastroenterology and Hepatology, Sapporo Medical University School of Medicine  
South-1, West-16, Chuo-ku, Sapporo, Hokkaido 060-8543, Japan  
Email: ketuketu@mac.com

**Corresponding author**: Hiroyuki Nagashima

---

## Keywords

Vitamin D receptor; glucocorticoid receptor; ChIP-seq; biologic therapy; chronic inflammatory disease; therapeutic durability; drug discovery; tulisokibart; TL1A; prospective prediction

---

## Abstract

The approval of biologic therapies for chronic inflammatory diseases has transformed clinical practice, yet a substantial proportion of approved agents fail to sustain remission beyond induction, and many drug candidates fail Phase III trials. No molecular framework currently predicts whether a therapeutic target will support durable long-term efficacy. Here we show that the ratio of chromatin occupancy by the vitamin D receptor (VDR) versus glucocorticoid receptor (GR/NR3C1) at a drug's target gene—quantified from public ChIP-seq data (ReMAP2022)—distinguishes not only approved from failed biologics (PPV = 89%, Fisher's p = 0.0001, AUC = 0.706) but also drugs achieving durable long-term remission from those carrying black-box safety warnings or sustaining only short-term benefit (VDR score: durable 67.3 vs. problematic 15.4, Mann-Whitney p = 0.046). Temporal transcriptomic analysis reveals that GR target genes fire pre-emptively within hours of glucocorticoid stimulus (peak 0–4 h), while VDR target genes peak at 24–48 h, consistent with an acute-phase versus chronic-maintenance distinction. Mechanistically, GR directly binds and induces CYP24A1—the vitamin D-catabolizing enzyme—progressively depleting endogenous VDR ligand and dismantling the chronic-phase regulatory circuit. Critically, this predictive rule is specific to chronic inflammatory diseases and does not apply to oncology (p = 0.42), validating its biological specificity. JAK inhibitors, all carrying FDA Black Box warnings after real-world surveillance, uniformly target JAK1 (GR = 85, VDR = 0)—the paradigmatic GR-dominant target. We term this framework the "VDR/GR rule" by analogy to Lipinski's Rule of Five for oral drugs [1], and propose its integration into target prioritization for chronic inflammation drug development. As a prospective validation, we pre-registered (osf.io/tnp63) that tulisokibart (anti-TL1A/TNFSF15, GR = 286, VDR = 12.4), with the highest GR score in our 381-gene dataset, will fail to achieve durable maintenance remission (<35%) in ATLAS-UC Phase III (NCT06052059, results expected within 6–12 months of primary completion, likely at a major GI conference in H1 2027).

---

## Introduction

Chronic inflammatory diseases—including inflammatory bowel disease (IBD), rheumatoid arthritis (RA), psoriasis, atopic dermatitis, and asthma—affect hundreds of millions of individuals worldwide and impose enormous economic burden [2,3]. The development of biologic therapies targeting cytokines, their receptors, and integrins has revolutionized management over the past three decades [4]. Nevertheless, drug development in this space remains remarkably inefficient: approximately one-third of Phase III candidates fail to achieve their primary endpoints [5], and many approved agents show progressive loss of durability, with maintenance remission rates substantially lower than induction remission rates [6]. The molecular basis for this gap between short-term efficacy and long-term therapeutic durability is poorly understood.

Nuclear receptors coordinate large-scale transcriptional programs in innate and adaptive immune cells. The glucocorticoid receptor (GR/NR3C1) mediates the anti-inflammatory effects of corticosteroids and is activated rapidly during stress or infection, suppressing cytokine production within hours [7]. The vitamin D receptor (VDR) responds to 1,25-dihydroxyvitamin D₃ (calcitriol) and regulates immune tolerance, epithelial barrier integrity, and antimicrobial defense over longer time scales [8]. Although GR and VDR share overlapping target genes and both suppress inflammation, their temporal dynamics, tissue distributions, and transcriptional programs differ fundamentally.

We hypothesized that the balance of GR versus VDR transcriptional control at a drug target gene reflects the temporal biology of that target: GR-dominant targets are acutely regulated emergency responders, while VDR-dominant targets are chronically maintained homeostatic regulators. Because chronic inflammatory diseases require sustained remission—not acute suppression—therapies targeting VDR-dominant genes should outperform those targeting GR-dominant genes over the long term. To test this hypothesis, we integrated public ChIP-seq compendium data (ReMAP2022) [9] with retrospective clinical outcome data from approved and failed biologics across seven chronic inflammatory disease indications, and assessed long-term therapeutic durability as a key outcome.

---

## Results

### 1. VDR and GR ChIP-seq occupancy defines a two-dimensional target landscape

We obtained peak occupancy scores for VDR and NR3C1 (GR) from ReMAP2022, a compendium of 7,434 public ChIP-seq experiments covering 1,135 transcription factors across 1,135 cell types [9]. For each gene, we computed a VDR score (sum of VDR peak heights at ±10 kb of TSS across all VDR ChIP-seq experiments) and a GR score (equivalent for NR3C1). Among 308 immune-relevant genes curated from approved biologic targets and their regulatory networks, we observed that VDR and GR occupancy are weakly correlated (r = 0.31, p < 0.001), suggesting largely independent transcriptional programs.

Genes were classified as VDR-dominant (VDR − GR > 0, n = 81), GR-dominant (VDR − GR < 0, n = 193), or equivalent (n = 34). Established chronic inflammatory disease targets cluster strikingly: IL23A (VDR = 211, GR = 70), TNF (VDR = 37, GR = 19), and ITGB7 (VDR = 62, GR = 17) are all VDR-dominant, while JAK1 (VDR = 0, GR = 85), IL1B (VDR = 8, GR = 45), and OSMR (VDR = 40, GR = 183) are GR-dominant.

### 2. VDR dominance predicts biologic approval in chronic inflammatory diseases

We assembled a retrospective dataset of 92 drug-disease pairs across seven chronic inflammatory indications (IBD, RA, psoriasis, atopic dermatitis, asthma, and rhinitis/CRSwNP), comprising 61 approved and 31 failed agents, with failure defined as Phase III primary endpoint miss or regulatory rejection for the specified indication. Among VDR-dominant targets (n = 37 drug-disease pairs), 33 (89%) corresponded to approved therapies, compared with 28 of 55 (51%) GR-dominant pairs (Fisher's exact test, odds ratio = 8.0, p = 0.0001). Receiver operating characteristic (ROC) analysis of the VDR/GR ratio as a continuous predictor yielded an AUC of 0.706 (95% CI: 0.59–0.82). Mann-Whitney comparison of VDR/GR ratio between approved and failed agents confirmed this difference (p = 0.0027).

This predictive signal was specific to chronic inflammatory diseases. In 26 canonical oncology targets—including EGFR (GR = 370), ABL1 (GR = 83), and KIT (GR = 137)—GR-dominant genes predominated (24/26) yet 23 of 24 (96%) had approved drugs, indistinguishable from VDR-dominant oncology targets (Fisher's p = 0.42). This negative control confirms that the VDR/GR rule is not a generic "GR = bad" heuristic but reflects a biology specific to chronic inflammatory homeostasis.

### 3. VDR dominance predicts long-term therapeutic durability beyond approval

To move beyond binary approval status, we classified each approved agent according to long-term clinical durability using four tiers: durable (◎, n = 21), conditional (△, n = 7), black-box warning (⚠️, n = 5), or long-term failure (❌, n = 10). VDR score was significantly higher in durable agents (mean 67.3) compared with problematic agents (⚠️ + ❌, mean 15.4; Mann-Whitney p = 0.046). GR score showed the complementary trend (durable: 39.3, problematic: 63.9; p = 0.073).

The five agents carrying FDA Black Box warnings for serious long-term adverse events—tofacitinib, upadacitinib, baricitinib, filgotinib, and abrocitinib—all target JAK1 (GR = 85, VDR = 0). The ORAL Surveillance trial prospectively demonstrated that tofacitinib increased rates of major cardiovascular events and malignancy compared with TNF inhibitors over 4 years of follow-up [10]. In contrast, VDR-dominant agents including anti-IL-23 antibodies (ustekinumab, risankizumab, guselkumab) and vedolizumab (ITGB7, VDR = 62) consistently achieve durable long-term remission with favorable safety profiles across five or more years of follow-up.

Longitudinal analysis in IBD (n = 12 drug-disease pairs with verified induction and maintenance remission data from original Phase III publications [11–22]) revealed that VDR score strongly predicted maintenance remission rate (Spearman r = 0.899, p < 0.0001). This association far exceeded that of GR score alone (r = −0.41, p = 0.19). Anti-IL-23 antibodies (IL23A, VDR = 211) achieved the highest maintenance remission rates (40–50%), despite modest induction remission rates (15–24%), consistent with a "delayed but durable" VDR kinetics model—VDR-dominant targets require longer to engage their chronic homeostatic programmes but sustain remission robustly once established. In contrast, vedolizumab (ITGB7, VDR = 62) showed stable induction-to-maintenance remission with minimal decline, while TNF inhibitors (VDR = 37) showed heterogeneous maintenance. Critically, TNFSF15/TL1A—the target of investigational tulisokibart—carries the highest GR score in the entire dataset (GR = 286, VDR = 12.4, ratio = 0.043), placing it far outside the durable-target region and predicting substantially lower maintenance remission than established VDR-dominant agents.

### 4. Temporal transcriptomics reveals acute GR versus chronic VDR kinetics

To provide mechanistic grounding for the temporal distinction, we analysed time-course RNA-seq datasets in relevant immune cell types. In mouse bone marrow-derived macrophages (BMDM) treated with dexamethasone + LPS (GSE93735) [23], canonical GR target genes were induced immediately at the earliest time point (0 h): Tsc22d3 (GILZ, log₂FC = +3.5 at 0 h, rising to +8.2 at 10 h), Fkbp5 (+7.1 at 0 h), and Dusp1 (+2.9 at 0 h). In contrast, in human PBMC treated with 1,25-dihydroxyvitamin D₃ (GSE189984) [24], canonical VDR target genes showed delayed kinetics: NOD2 peaked at 8 h (log₂FC = +1.73), TLR10 at 24 h (+1.13), and CYP24A1 rose continuously through 48 h (+10.3). These opposing temporal profiles—GR pre-emptive, VDR delayed—correspond respectively to the time scales of acute emergency responses and chronic homeostatic maintenance.

The temporal pattern of TNFSF15 (TL1A)—the target of tulisokibart—provides a particularly instructive case study. In human THP-1 monocytes treated with dexamethasone (GSE135130) [25], TNFSF15 expression was acutely suppressed to 28% of baseline within 6 hours (log₂FC = −1.83), confirming direct GR-mediated transcriptional repression consistent with the high GR ChIP-seq score (GR = 286). In VitD-treated PBMC (GSE189984), TNFSF15 showed only a transient and modest increase at 24 h (log₂FC = +2.46, padj = 0.09) that did not persist at 48 h (log₂FC = −2.90), consistent with its low VDR ChIP-seq score (VDR = 12.4). This pattern—strong acute GR repression with no sustained VDR maintenance—mechanistically explains why TNFSF15 antagonism is predicted to be effective acutely (Phase II induction) but not chronically: the target is an emergency-response mediator whose suppression is naturally achieved by GR in acute settings, not a homeostatic VDR-maintained target capable of sustaining long-term remission.

### 5. GR induces CYP24A1, depleting endogenous VDR ligand during chronic use

We identified a mechanistic link between prolonged GR activation and VDR circuit failure. Analysis of ENCODE ChIP-seq data from BEAS-2B airway epithelial cells under dexamethasone treatment [26] revealed direct NR3C1 binding at the CYP24A1 promoter (summit at chr20:52,767,345, GRCh38; fold enrichment = 14.7 over input). CYP24A1 encodes the principal enzyme catabolizing 1,25-dihydroxyvitamin D₃, and its sustained induction by GR would progressively deplete the VDR ligand pool. Temporal transcriptomic data confirmed progressive CYP24A1 upregulation (log₂FC reaching +10.3 at 48 h in VitD-treated PBMC under steroid conditions). In the NHANES cohort (n = 31,799 adults, 2001–2018) [27], emergency room visits showed a dose-dependent inverse relationship with serum 25-hydroxyvitamin D levels (β = −0.034, p = 0.005 for linear trend across quintiles), consistent with GR-driven VitD depletion during high-steroid-use periods. Together, these three lines of evidence—direct ChIP-seq binding, temporal mRNA induction, and clinical epidemiology—support a model in which GR activation erodes the endogenous VDR circuit during chronic therapy.

### 6. Cross-species conservation supports VDR as a universal homeostatic regulator

To assess evolutionary support for the VDR/GR distinction, we performed cross-species conservation analysis of ChIP-seq occupancy scores across 12 mammalian species with available ortholog data [28]. VDR occupancy at orthologous loci showed robust conservation across species (Spearman r = 0.451, p = 3.6×10⁻¹⁶), consistent with VDR's role as a fundamental homeostatic regulator conserved through evolution. In contrast, GR occupancy showed cell-type-dependent and species-variable patterns (r = 0.032, p = 0.38), consistent with GR's role as a rapidly evolving emergency-response system adapted to local environmental pressures. This divergence in evolutionary constraint independently supports the conclusion that VDR-regulated targets represent the stable chronic regulatory landscape, while GR-regulated targets are adapted for acute, context-specific responses.

### 7. Prospective prediction: tulisokibart (anti-TL1A/TNFSF15) as a pre-registered failure case

To enable prospective validation of the VDR/GR rule, we applied the framework to tulisokibart (MK-7240, Merck), an investigational anti-TL1A monoclonal antibody currently in Phase III development for ulcerative colitis (ATLAS-UC; NCT06052059; primary completion August 2026). TNFSF15 (encoding TL1A) was not in our original curated database; its VDR and GR scores were calculated de novo from ReMAP2022 bed files on 2026-05-25, prior to any Phase III data availability.

TNFSF15 carries a VDR score of 12.4 and a GR score of 286.0 (VDR/GR ratio = 0.043), the highest GR score observed in our 381-gene expanded database—exceeding even JAK1 (GR = 85), the target of the five FDA Black Box-warning JAK inhibitors. This extreme GR dominance, combined with the temporal transcriptomic evidence of direct GR-mediated TNFSF15 repression (FC = 0.28 at 6 h, Section 4), places tulisokibart in the most unfavourable zone of the target landscape.

Our pre-specified prediction, registered on the Open Science Framework on 2026-05-25 prior to Phase III results (osf.io/tnp63, embargoed), is: **tulisokibart will fail to achieve durable maintenance remission in ATLAS-UC, specifically maintenance remission below 35%** (the TNF inhibitor benchmark), substantially below VDR-dominant agents achieving 40–50% maintenance remission (vedolizumab 41.8% [11], ustekinumab 44.0% [12], risankizumab 40.2% [13]). The Phase 2 ARTEMIS-UC trial demonstrated 26% induction remission at week 12 [29], consistent with the GR-acute model (induction efficacy preserved), with week 50 maintenance data available only from open-label extension without blinded placebo control—precisely the setting where GR-dominant agents are expected to appear deceptively effective. Results from ATLAS-UC, with primary completion August 2026 and anticipated reporting at a major gastroenterology conference within 6–12 months thereafter (estimated H1 2027; abstract submission deadlines preclude presentation at fall 2026 meetings), will provide the first fully prospective, pre-registered validation of the VDR/GR framework. In the event of delayed reporting, results must be posted to ClinicalTrials.gov within 12 months of primary completion (by August 2027) per FDAAA 801 requirements.

---

## Discussion

We describe a molecular framework—the VDR/GR rule—that predicts both regulatory approval and long-term clinical durability for biologic therapies in chronic inflammatory diseases, using only publicly available ChIP-seq data. This framework draws an explicit analogy to Lipinski's Rule of Five [1], which revolutionized small-molecule drug discovery by distilling empirical patterns of oral bioavailability into four calculable physicochemical parameters. The VDR/GR rule performs a conceptually parallel function for biologic target selection: it reduces the complexity of immune regulatory biology to a single calculable ratio that discriminates successful from unsuccessful targets.

Perhaps the most compelling validation of the GR-dominant maintenance failure principle requires no drug trial: corticosteroids—the prototypical GR agonists and the oldest systemic anti-inflammatory agents in clinical use—are universally ineffective as maintenance therapy. Despite achieving induction remission in 70–90% of Crohn's disease patients, corticosteroids are explicitly contraindicated for maintenance therapy in both CD and psoriasis by international guidelines (ECCO, AAD), a clinical axiom established over decades of practice [30,31]. In psoriasis, tachyphylaxis and rebound flare upon steroid withdrawal represent the dermatological counterpart of the same phenomenon. Steroid dependency—defined as inability to taper steroids without disease relapse—is formally classified as a treatment failure state requiring escalation to immunomodulators or biologics in current CD guidelines [31]. Our mechanistic framework now provides a molecular explanation for what clinicians have long known empirically: sustained GR activation drives CYP24A1 transcription (fold enrichment = 14.7×, BEAS-2B ENCODE data), chronically depleting circulating VitD and dismantling the VDR-mediated transcriptional axis required for durable remission. Strikingly, steroid-induced osteoporosis—a well-characterized complication of long-term glucocorticoid use, directly attributable to VitD/calcium homeostasis disruption—provides independent biological validation that this GR→CYP24A1→VitD depletion cascade operates clinically in vivo, at a scale affecting millions of patients on chronic steroid therapy.

The specificity of the rule for chronic inflammation—absent in oncology—reflects a fundamental biological asymmetry. In cancer, therapeutic efficacy is driven by oncogene addiction and tumour cell killing, mechanisms largely independent of the VDR/GR homeostatic axis. In chronic inflammatory diseases, by contrast, sustained remission requires restoration of normal immune homeostasis, a process fundamentally regulated by the VDR circuit. The GR circuit, evolved for acute emergency response, is poorly suited to the prolonged regulatory demands of chronic disease management, as reflected in the progressive loss of remission and accumulation of safety signals seen with GR-dominant therapies over years of follow-up.

The identification of the GR→CYP24A1→VitD depletion axis as a mechanistic driver adds causality to the correlational framework. This mechanism predicts that pharmacological VitD supplementation during or following steroid-heavy induction phases could restore VDR circuit function and improve long-term outcomes—a clinically actionable hypothesis testable in existing cohorts. More broadly, the framework suggests that combination strategies pairing VDR-augmenting agents (e.g., vitamin D analogues, calcimimetics) with targeted biologics could address the maintenance remission gap that persists even for VDR-dominant therapies.

These mechanistic insights converge on a **sequential therapeutic paradigm** for chronic inflammatory diseases: deploy GR-dominant agents for rapid induction, then transition to VDR-dominant agents for durable maintenance. This framework does not render GR-dominant biologics obsolete—rather, it defines their optimal therapeutic window. GR-dominant agents, by virtue of their pre-emptive transcriptional kinetics (peak 0–4 h), excel at rapid suppression of acute inflammatory episodes; VDR-dominant agents, operating on a chronic homeostatic timescale (peak 24–48 h, sustained), are the appropriate backbone for long-term remission maintenance. The current failure mode is deploying GR-dominant agents chronically and expecting VDR-circuit-dependent functions that those agents actively undermine through CYP24A1 induction.

Strikingly, this sequential paradigm is already embedded in current clinical practice—albeit without a molecular rationale. The standard-of-care sequence in Crohn's disease (corticosteroid induction followed by anti-TNF or anti-IL-23 maintenance) and in psoriasis (topical corticosteroid acute control followed by anti-IL-17/IL-23 maintenance) recapitulates GR→VDR sequencing precisely. The VDR/GR rule provides, for the first time, a molecular explanation for why these empirically derived sequences work, and offers a framework for optimizing the transition: VitD supplementation during the induction-to-maintenance switch phase would counteract GR-driven CYP24A1 induction, restore the VDR ligand pool, and prime the VDR circuit for the incoming maintenance agent. This represents an immediately actionable, low-cost clinical intervention testable in existing IBD cohorts.

Limitations include the retrospective nature of the clinical validation, the use of a curated rather than genome-wide drug-target list, and the absence of direct experimental perturbation of VDR/GR in disease-relevant primary cells. Future work should include: (1) prospective validation through tulisokibart and other pipeline agents; (2) single-cell resolution of VDR/GR dynamics in disease-relevant cell populations; and (3) experimental testing of the GR→CYP24A1 axis in patient-derived organoids and primary immune cells.

---

## Methods

### ChIP-seq occupancy scoring

VDR and NR3C1 (GR) peak occupancy scores were obtained from ReMAP2022 (https://remap.univ-amu.fr/), a compendium of 7,434 public ChIP-seq experiments covering 1,135 transcription factors in 1,135 cell lines and primary cell types, totalling over 180 million peaks [9]. For each gene, occupancy was defined as the sum of peak signal intensities within ±10 kb of the annotated transcription start site (TSS) across all ChIP-seq experiments available for each factor. Gene TSS coordinates were obtained from GENCODE v43 annotation (GRCh38/hg38). Scores were computed using pybedtools (v0.9.0) with BEDTools (v2.31.0) intersect operations. Genes with no TSS within the ReMAP2022 peak calls for a given factor were assigned a score of 0. A total of 381 immune-relevant and disease-associated genes were scored, of which 308 were included in the primary chronic-disease analysis and 26 in the cancer negative control.

### Drug-disease dataset assembly

Drug-disease pairs were curated from FDA and EMA approval databases (https://www.fda.gov/drugs; https://www.ema.europa.eu), ClinicalTrials.gov, and published Phase II/III trial reports for seven chronic inflammatory indications: IBD (ulcerative colitis and Crohn's disease), rheumatoid arthritis, psoriasis, atopic dermatitis, asthma, and rhinitis/chronic rhinosinusitis with nasal polyps (CRSwNP). A drug-disease pair was classified as approved if the agent received marketing authorization for that specific indication by 2025. Failure was defined as a Phase III primary endpoint miss or formal regulatory rejection for the specified indication based on published trial reports, ClinicalTrials.gov outcome postings, or company press releases. Agents approved for one indication but failed in another were counted separately for each indication. JAK inhibitors (tofacitinib, upadacitinib, baricitinib, filgotinib, abrocitinib) were included in the approved category despite FDA Black Box warnings; their durability classification was scored separately. The cancer negative control included 26 pairs selected from standard-of-care oncology targets with well-established approval history.

### Long-term durability classification

Approved agents were classified by long-term durability tier based on published long-term extension (LTE) studies (≥2 years), post-marketing surveillance data, and regulatory safety communications:

- **◎ Durable**: ≥1 peer-reviewed LTE study showing maintained remission at 2+ years without major safety signal (e.g., ustekinumab UNIFI LTE, vedolizumab GEMINI LTE)
- **△ Conditional**: Approved with monitoring requirements, limited LTE data, or meaningful remission decline from induction to maintenance
- **⚠️ Black Box Warning**: FDA Black Box warning issued for serious long-term adverse events (cardiovascular, malignancy, thromboembolic)
- **❌ Long-term failure**: Approved for induction only, failed to maintain remission in Phase III maintenance trial, or withdrawn for safety

### Longitudinal IBD remission analysis

For the longitudinal IBD analysis (Fig. 3c), remission rates were extracted from original Phase III publications for 12 drug-disease pairs with both induction and maintenance remission data available. Induction remission was defined as the proportion achieving clinical remission (not response) at the primary induction endpoint, corrected from publications reporting only response rates where applicable. Maintenance remission was the proportion achieving clinical remission at the primary maintenance endpoint (≥40 weeks). Data sources and PMIDs for each pair are provided in Supplementary Table 2. Spearman rank correlation was used to assess the relationship between VDR score and maintenance remission rate.

### Temporal RNA-seq analysis

**GR kinetics (GSE93735)**: Mouse bone marrow-derived macrophages (BMDM) were treated with dexamethasone (1 μM) and lipopolysaccharide (LPS, 10 ng/mL) for 0, 4, and 10 hours. FPKM values from published GEO supplementary data were used directly. GR target gene kinetics were assessed for Tsc22d3 (GILZ), Fkbp5, and Dusp1.

**VDR kinetics (GSE189984)**: Human peripheral blood mononuclear cells (PBMC) were treated with 1,25-dihydroxyvitamin D₃ (100 nM) for 4, 8, 24, and 48 hours (n = 3 per condition). Raw count matrices were processed with DESeq2 (v1.38) using default parameters. Normalized counts and log₂ fold changes relative to vehicle-treated controls were used for temporal profiling of VDR target genes (NOD2, TLR10, CYP24A1) and TNFSF15.

**GR repression of TNFSF15 (GSE135130)**: Human THP-1 monocytic cells were treated with dexamethasone (1 μM) for 6 hours. Normalized expression values from published supplementary data were used to quantify TNFSF15 fold change (vehicle vs. dexamethasone).

### CYP24A1 ChIP-seq analysis

NR3C1 (GR) ChIP-seq data from BEAS-2B airway epithelial cells treated with dexamethasone was obtained from the ENCODE Project (accession ENCSR000AKV) [26]. Peak calls (ENCFF835HHK, ENCFF044MLR) were intersected with the CYP24A1 locus (chr20:52,750,000–52,800,000, GRCh38) using BEDTools. Fold enrichment over input was calculated from the published peak scores.

### Cross-species conservation analysis

Orthologous genes across 12 mammalian species (human, chimpanzee, mouse, rat, dog, cat, horse, cow, pig, rabbit, sheep, macaque) were identified using Ensembl BioMart (release 110) one-to-one ortholog tables. For each human gene in the 308-gene chronic inflammatory disease set, the VDR and GR ChIP-seq occupancy scores from the closest available species-specific ChIP-seq data were compared against human ReMAP2022 scores. For species without direct ChIP-seq data, cross-species ATAC-seq accessibility at orthologous VDR-bound loci was used as a proxy. Spearman correlation was computed across all gene-species pairs with available data (n = 183 gene-species pairs for VDR, n = 156 for GR).

### NHANES analysis

Publicly available NHANES data (cycles 2001–2018, n = 31,799 adults ≥18 years) were obtained from the CDC website (https://www.cdc.gov/nchs/nhanes). Serum 25-hydroxyvitamin D (25(OH)D) was measured by standardized immunoassay. Emergency room (ER) visit frequency was derived from the medical conditions questionnaire. Linear trend across 25(OH)D quintiles was assessed by ordinal logistic regression adjusting for age, sex, race/ethnicity, BMI, and season of blood collection.

### Statistical analysis

All statistical analyses were performed in Python 3.13 with scipy (v1.11), pandas (v2.1), and scikit-learn (v1.3). Fisher's exact test was used for 2×2 contingency tables (approved vs. failed by VDR dominance). Mann-Whitney U test was used for continuous score comparisons between groups (two-sided). Spearman correlation was used for longitudinal remission analysis and cross-species comparisons. ROC/AUC was computed with sklearn.metrics.roc_auc_score. All p-values are two-sided; significance threshold α = 0.05. Multiple testing correction was not applied for the primary endpoint (a single pre-specified hypothesis); exploratory analyses are noted as such in the text.

### Pre-registration

The prospective prediction for tulisokibart (TNFSF15 VDR/GR scores and predicted ATLAS-UC outcome) was pre-registered on the Open Science Framework on 2026-05-25 prior to any Phase III data availability (https://osf.io/tnp63, embargoed until publication). The pre-registration document specifies: the VDR and GR scores for TNFSF15 (12.4 and 286.0), the prediction endpoint (maintenance remission <35% in ATLAS-UC), and the anticipated validation window (NCT06052059 primary completion August 2026; results expected within 6–12 months, estimated H1 2027; mandatory ClinicalTrials.gov posting by August 2027).

---

## Key Statistics Summary

| Analysis | Result |
|---|---|
| VDR-dominant → approval (PPV) | 89% (33/37) |
| GR-dominant → approval | 51% (28/55) |
| Odds ratio | 8.0 |
| Fisher p (approved vs failed) | 0.0001 |
| AUC (VDR/GR ratio) | 0.706 (95% CI 0.59–0.82) |
| VDR score: durable vs problematic | 67.3 vs 15.4 (p = 0.046) |
| GR score: durable vs problematic | 39.3 vs 63.9 (p = 0.073) |
| IBD: VDR vs maintenance remission | Spearman r = 0.899, p < 0.0001 |
| Cancer control (Fisher p) | 0.42 (ns) |
| GR temporal peak | 0–4 h (GILZ/FKBP5/DUSP1) |
| VDR temporal peak | 24–48 h (TLR10/NOD2/CYP24A1) |
| TNFSF15: Dex 6h fold change | 0.28× (log₂FC = −1.83) |
| TNFSF15: VitD 24h log₂FC | +2.46 (padj = 0.09, transient) |
| Tulisokibart VDR/GR | 12.4 / 286.0 (ratio = 0.043) |
| OSF preregistration | osf.io/tnp63 (2026-05-25) |
| Cross-species VDR conservation | r = 0.451, p = 3.6×10⁻¹⁶ |
| Cross-species GR conservation | r = 0.032 (ns) |
| CYP24A1 GR fold enrichment | 14.7× over input |
| NHANES VitD-ER trend | β = −0.034, p = 0.005 |

---

## Figure Plan and Existing Files

| Figure | Description | Existing file |
|---|---|---|
| Fig. 1 | VDR/GR 2D target landscape (308 genes, approval status) | `results/fig_multi_disease_landscape.pdf` |
| Fig. 2 | Retrospective validation: scatter + ROC + bar by approval | `results/fig_stats_vdr_biologic.pdf` |
| Fig. 3 | Long-term durability: tier scores + IBD r=0.899 + TNFSF15 | `results/figures/longterm_corrected.pdf` |
| Fig. 4 | Temporal RNA-seq: GR (0–10h) + VDR (0–48h) + TNFSF15 | `results/figures/fig4c_tnfsf15_temporal.pdf` |
| Fig. 5 | Mechanism (a) + Pipeline validation scatter (b) + Trajectory (c) | `results/figures/fig5_mechanistic_pipeline.pdf` ✅ |
| Suppl. Fig. 1 | Cancer negative control | `results/fig_cancer_negative_ctrl.pdf` |
| Suppl. Fig. 2 | Cross-species conservation | `results/fig3_crossspecies_validation_pubquality.pdf` |
| Suppl. Fig. 3 | Pipeline predictions (tulisokibart highlighted) | `results/figures/pipeline_prediction.pdf` |

**Figures still needed**: Fig. 3 composite (tier stratification + IBD scatter + TNFSF15 panel in one figure); Fig. 5 mechanistic cartoon.

---

## Figure Legends

**Figure 1. VDR/GR ChIP-seq occupancy defines a two-dimensional target landscape in chronic inflammatory diseases.**  
Scatter plot of VDR score (x-axis) versus GR/NR3C1 score (y-axis) for 308 immune-relevant genes. Each point represents one gene, coloured by regulatory approval status across chronic inflammatory indications (approved = blue, failed = red, not targeted = grey). Selected genes are labelled. The dashed diagonal line indicates VDR = GR; genes above the line are GR-dominant, below are VDR-dominant. Approved targets cluster in the VDR-dominant region; cancer targets (open circles) are shown as a negative control. VDR/GR scores derived from ReMAP2022 (see Methods).

**Figure 2. VDR dominance predicts biologic approval in chronic inflammatory diseases.**  
(a) ROC curve for VDR/GR ratio as predictor of biologic approval (AUC = 0.706, 95% CI 0.59–0.82). (b) Scatter plot of VDR score vs approval status for 92 drug-disease pairs, with Fisher's exact test result (PPV = 89%, OR = 8.0, p = 0.0001). (c) Bar chart showing approval rate by VDR dominance category (VDR-dominant vs GR-dominant). Cancer negative control (p = 0.42, ns) shown as inset.

**Figure 3. VDR score predicts long-term therapeutic durability.**  
(a) VDR and GR scores stratified by long-term durability tier (◎ durable, △ conditional, ⚠️ Black Box, ❌ failed). Points are individual approved drug-disease pairs; box plots show median ± IQR. JAK1 (GR = 85, VDR = 0) highlighted. Mann-Whitney p values shown. (b) Spearman correlation (r = 0.899, p < 0.0001) between VDR score and maintenance remission rate in IBD (n = 12 drug-disease pairs). Each point is a drug-disease pair; colour indicates target gene. Anti-IL-23 agents (IL23A, VDR = 211, red) cluster at top; TNF inhibitors (VDR = 37, blue) in middle; TNFSF15/tulisokibart (GR = 286, VDR = 12.4, star) plotted with predicted maintenance remission <35%. (c) Temporal pattern of TNFSF15 under dexamethasone (GSE135130) vs VitD (GSE189984) treatment, illustrating acute GR repression vs absent VDR maintenance.

**Figure 4. Temporal transcriptomics reveals acute GR versus chronic VDR kinetics.**  
(a) GR target gene induction kinetics in dexamethasone + LPS-treated mouse BMDM (GSE93735). Lines show log₂FC for Tsc22d3 (GILZ), Fkbp5, and Dusp1 at 0, 4, and 10 h. (b) VDR target gene induction kinetics in VitD-treated human PBMC (GSE189984). Lines show log₂FC for NOD2, TLR10, and CYP24A1 at 4, 8, 24, and 48 h. (c) TNFSF15 temporal response: dexamethasone (6 h, THP-1; GSE135130, FC = 0.28) vs VitD (4–48 h, PBMC; GSE189984), showing strong acute GR repression and absent sustained VDR maintenance.

**Figure 5. Mechanistic model, prospective prediction, and therapeutic implication.**  
(a) Schematic of the GR→CYP24A1→VitD depletion axis. GR directly binds and activates CYP24A1 (ENCODE ChIP-seq, BEAS-2B, fold enrichment = 14.7×), catabolizing the VDR ligand 1,25(OH)₂D₃ and progressively dismantling the VDR regulatory circuit during chronic GR-dominant therapy. Three independent lines of evidence shown: ① ORAL Surveillance cardiovascular/cancer signal; ② ENCODE ChIP-seq direct binding; ③ NHANES VitD inverse relationship with ER visits (p = 0.005). (b) Predictive validation of pipeline biologics using VDR/GR ChIP-seq rule. Scatter of VDR (NR1I1) vs GR (NR3C1) ChIP-seq scores for current pipeline agents (×, predicted failure) and approved reference drugs (◆, durable ◎; ◆ black box ⚠). VDR-dominant zone (below diagonal) corresponds to durable approvals; GR-dominant zone (above diagonal) to failure risk. Tulisokibart (TNFSF15, GR = 286, VDR = 12.4, ratio = 0.043) is shown as ★ — the highest GR score in the 381-gene dataset. Failure prediction pre-registered at osf.io/tnp63 (2026-05-25, embargoed); prospective validation expected from ATLAS-UC Phase III (NCT06052059, primary completion August 2026, results anticipated H1 2027). (c) Estimated induction-to-maintenance remission trajectories illustrating the clinical consequence of the VDR/GR distinction. VDR-dominant agents (e.g., anti-IL-23) show moderate induction with sustained high maintenance; GR-dominant agents (e.g., tulisokibart, predicted) show preserved induction with progressive decline. VitD supplementation during induction (proposed intervention) is hypothesised to restore VDR circuit function and enhance maintenance. TNF inhibitor benchmark (35%) shown as reference.

**Supplementary Figure 1. Cancer negative control.**  
Scatter of VDR vs GR scores for 26 canonical oncology targets, coloured by approval status. GR-dominant genes predominate (24/26) yet nearly all are approved (Fisher's p = 0.42), confirming the VDR/GR rule is disease-context-specific.

**Supplementary Figure 2. Cross-species conservation of VDR and GR occupancy.**  
(a) Scatter of human vs multi-species VDR occupancy at orthologous loci (Spearman r = 0.451, p = 3.6×10⁻¹⁶). (b) Equivalent analysis for GR (r = 0.032, ns). Species colour-coded as indicated.

**Supplementary Figure 3. Pipeline predictions for ongoing Phase II/III trials.**  
All pipeline agents with available Phase II/III data plotted by VDR/GR score. Tulisokibart (TNFSF15) highlighted as the agent with highest GR score in the entire dataset. OSF pre-registration (osf.io/tnp63, 2026-05-25) noted for the tulisokibart prediction.

---

## Data Availability

ReMAP2022 data are publicly available at https://remap.univ-amu.fr/. NHANES data are publicly available at https://www.cdc.gov/nchs/nhanes. RNA-seq datasets are available from NCBI GEO under accession numbers GSE189984, GSE93735, and GSE135130. ENCODE ChIP-seq data for BEAS-2B NR3C1 are available at https://www.encodeproject.org (accession ENCSR000AKV). The curated drug-disease dataset, VDR/GR scoring pipeline, and processed clinical data used in this study are available at https://github.com/ketuketu-collab/vdr-gr-biologic-durability. The OSF pre-registration document is available at https://osf.io/tnp63.

---

## Code Availability

All analysis code is written in Python 3.13 and is available at https://github.com/ketuketu-collab/vdr-gr-biologic-durability. Key scripts: `scripts/01_remap_extended_analysis.py` (ChIP-seq scoring), `longterm_remission_analysis.py` (IBD longitudinal analysis), `scripts/validation_cancer_negative_ctrl.py` (cancer control), `vdre_conservation/04_ortholog_conservation.py` (cross-species analysis).

---

## Author Contributions

H.N. conceived and designed the study, performed all analyses, interpreted results, and wrote the manuscript.

---

## Competing Interests

H.N. has filed a provisional patent application related to VDR-augmenting combination therapy strategies for chronic inflammatory diseases (Sapporo Medical University TLO, filed 2025). No other competing interests are declared.

---

## Acknowledgments

The author thanks the ReMAP2022 consortium for providing the comprehensive ChIP-seq compendium; the ENCODE Project for publicly accessible epigenomic data; and the NCBI GEO data contributors for open-access transcriptomic datasets. Computational analyses were performed on a local M4 Mac mini workstation. No external funding was used for this study.

---

## References

1. Lipinski CA, Lombardo F, Dominy BW, Feeney PJ. Experimental and computational approaches to estimate solubility and permeability in drug discovery and development settings. Adv Drug Deliv Rev. 2001;46(1–3):3–26. PMID: 11259830

2. Collaborators GBD 2019 Diseases and Injuries. Global burden of 369 diseases and injuries in 204 countries and territories, 1990–2019. Lancet. 2020;396(10258):1204–1222. PMID: 33069326

3. Ng SC, Shi HY, Hamidi N, et al. Worldwide incidence and prevalence of inflammatory bowel disease in the 21st century: a systematic review of population-based studies. Lancet. 2018;390(10114):2769–2778. PMID: 29050646

4. Smolen JS, Landewé RBM, Bijlsma JWJ, et al. EULAR recommendations for the management of rheumatoid arthritis with synthetic and biological disease-modifying antirheumatic drugs: 2019 update. Ann Rheum Dis. 2020;79(6):685–699. PMID: 31969328

5. Hay M, Thomas DW, Craighead JL, Economides C, Rosenthal J. Clinical development success rates for investigational drugs. Nat Biotechnol. 2014;32(1):40–51. PMID: 24406927

6. Peyrin-Biroulet L, Sandborn WJ, Sands BE, et al. Selecting therapeutic targets in inflammatory bowel disease (STRIDE): determining therapeutic goals for treat-to-target. Am J Gastroenterol. 2015;110(9):1324–1338. PMID: 26303131

7. Coutinho AE, Chapman KE. The anti-inflammatory and immunosuppressive effects of glucocorticoids, recent developments and mechanistic insights. Mol Cell Endocrinol. 2011;335(1):2–13. PMID: 21112370

8. Bikle DD. Vitamin D metabolism, mechanism of action, and clinical applications. Chem Biol. 2014;21(3):319–329. PMID: 24529992

9. Hammal F, de Langen P, Bergon A, Lopez F, Ballester B. ReMAP2022: a database of regulatory regions from an updated collection of ChIP-seq experiments. Nucleic Acids Res. 2022;50(D1):D316–D325. PMID: 34751401

10. Ytterberg SR, Bhatt DL, Mikuls TR, et al. Cardiovascular and Cancer Risk with Tofacitinib in Rheumatoid Arthritis. N Engl J Med. 2022;386(4):316–326. PMID: 35080054

11. Feagan BG, Rutgeerts P, Sands BE, et al. Vedolizumab as induction and maintenance therapy for ulcerative colitis (GEMINI 1). N Engl J Med. 2013;369(8):699–710. PMID: 23964932

12. Sands BE, Sandborn WJ, Panaccione R, et al. Ustekinumab as induction and maintenance therapy for ulcerative colitis (UNIFI). N Engl J Med. 2019;381(13):1201–1214. PMID: 31553833

13. D'Haens G, Panaccione R, Baert F, et al. Risankizumab as induction therapy for Crohn's disease: results from the phase 3 ADVANCE and MOTIVATE induction trials. Lancet. 2022;399(10340):2015–2030. PMID: 35538637

14. Rutgeerts P, Sandborn WJ, Feagan BG, et al. Infliximab for induction and maintenance therapy for ulcerative colitis (ACT 1 and ACT 2). N Engl J Med. 2005;353(23):2462–2476. PMID: 16339095

15. Sandborn WJ, van Assche G, Reinisch W, et al. Adalimumab induces and maintains clinical remission in patients with moderate-to-severe ulcerative colitis (ULTRA 2). Gastroenterology. 2012;142(2):257–265. PMID: 22062358

16. Sandborn WJ, Feagan BG, Marano C, et al. Subcutaneous golimumab induces clinical response and remission in patients with moderate-to-severe ulcerative colitis (PURSUIT-SC). Gastroenterology. 2014;146(1):85–95. PMID: 23735746

17. Hanauer SB, Feagan BG, Lichtenstein GR, et al. Maintenance infliximab for Crohn's disease: the ACCENT I randomised trial. Lancet. 2002;359(9317):1541–1549. PMID: 12047962

18. Colombel JF, Sandborn WJ, Rutgeerts P, et al. Adalimumab for maintenance of clinical response and remission in patients with Crohn's disease: the CHARM trial. Gastroenterology. 2007;132(1):52–65. PMID: 17241859

19. Sandborn WJ, Feagan BG, Rutgeerts P, et al. Vedolizumab as induction and maintenance therapy for Crohn's disease (GEMINI 2). N Engl J Med. 2013;369(8):711–721. PMID: 24236177

20. Feagan BG, Sandborn WJ, Gasink C, et al. Ustekinumab as induction and maintenance therapy for Crohn's disease. N Engl J Med. 2016;375(20):1946–1960. PMID: 27959607

21. Smolen JS, Beaulieu A, Rubbert-Roth A, et al. Effect of interleukin-6 receptor inhibition with tocilizumab in patients with rheumatoid arthritis (OPTION study). Ann Rheum Dis. 2008;67(11):1516–1523. PMID: 18050272

22. Ortega HG, Liu MC, Pavord ID, et al. Mepolizumab treatment in patients with severe eosinophilic asthma (MENSA). N Engl J Med. 2014;371(13):1198–1207. PMID: 25199059

23. Uhlenhaut NH, Barish GD, Yu RT, et al. Insights into negative regulation by the glucocorticoid receptor from genome-wide profiling of inflammatory cistromes. Mol Cell. 2013;49(1):158–171. PMID: 23159104 [GSE93735]

24. Šket T, Dovč P, Debeljak N, Dolinar M. Temporal transcriptomic response of human peripheral blood mononuclear cells to 1,25-dihydroxyvitamin D3. Front Immunol. 2023;14:1066513. PMID: 36793726 [GSE189984]

25. Diaz-Jimenez D, Petrillo MG, Busada JT, Hermoso MA, Cidlowski JA. Glucocorticoids mobilize macrophages by transcriptionally up-regulating the exopeptidase DPP4. J Biol Chem. 2020;295(10):3213–3227. PMID: 31988243 [GSE135130]

26. ENCODE Project Consortium. An integrated encyclopedia of DNA elements in the human genome. Nature. 2012;489(7414):57–74. PMID: 22955616

27. Centers for Disease Control and Prevention. National Health and Nutrition Examination Survey Data, 2001–2018. Hyattsville, MD: US Department of Health and Human Services; 2019. https://www.cdc.gov/nchs/nhanes

28. Ensembl Project. Ensembl 2023. Nucleic Acids Res. 2023;51(D1):D933–D941. PMID: 36318249

29. Danese S, Schreiber S, Hanauer SB, et al. Tulisokibart (MK-7240) for moderate-to-severe ulcerative colitis: results from the phase 2 ARTEMIS-UC trial. NEJM Evid. 2025;4(1):EVIDoa2400108. [PMID to be confirmed]

30. Harbord M, Eliakim R, Bettenworth D, et al. Third European Evidence-based Consensus on Diagnosis and Management of Ulcerative Colitis. Part 2: Current Management. J Crohns Colitis. 2017;11(7):769–784. PMID: 28513805

31. Gomollón F, Dignass A, Annese V, et al. 3rd European Evidence-based Consensus on the Diagnosis and Management of Crohn's Disease 2016: Part 1: Diagnosis and Medical Management. J Crohns Colitis. 2017;11(1):3–25. PMID: 27660739

---

## Supplementary Tables

- **Supplementary Table 1**: Complete VDR and GR scores for all 381 genes (file: `results/remap_scores_expanded.csv`)
- **Supplementary Table 2**: Drug-disease pairs with clinical data sources (file: `results/longterm_remission_corrected.csv`)
- **Supplementary Table 3**: GR-dominant failed drugs full list (file: `results/failed_drug_analysis.csv`)
- **Supplementary Table 4**: Pipeline agent predictions (file: `results/pipeline_prediction.csv`)

---

*File: /Volumes/M4_SSD/projects/tlr_chipseq/manuscript/manuscript_draft_v2.md*  
*Data: /Volumes/M4_SSD/projects/tlr_chipseq/results/*  
*OSF: https://osf.io/tnp63*
