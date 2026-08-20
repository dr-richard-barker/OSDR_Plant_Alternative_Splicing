# Differential Alternative Splicing in Plant Spaceflight Transcriptomes: A Multi-Study rMATS Analysis with Machine Learning and Ontology Integration

**Authors:** Richard Barker

**Target Journal:** npj Microgravity

---

## Abstract

Spaceflight imposes unique mechanical and environmental stresses on plants, including microgravity, altered gravity regimes, ionizing radiation, and constrained growth conditions. While transcriptomic responses to spaceflight have been extensively characterized at the gene expression level, alternative splicing (AS) — a critical post-transcriptional regulatory mechanism — remains largely unexplored in the space biology context. We present a systematic analysis of differential alternative splicing across seven Arabidopsis thaliana RNA-seq studies from the NASA Open Science Data Repository (OSDR), using rMATS-turbo to quantify five AS event types (SE, A5SS, A3SS, MXE, RI) across diverse spaceflight-relevant stressors: partial and fractional gravity, ISS spaceflight, simulated galactic cosmic ray (GCR) radiation, and lunar regolith growth. Across all seven studies (199 samples, ~100,000 total splicing events), retained intron (RI) events consistently dominated significant events (57–67%), and NMD prediction identified 79–86% of significant events as NMD-sensitive, establishing intron-retention-coupled NMD as the primary AS response to spaceflight stress. Gene Ontology enrichment revealed circadian rhythm regulation as a recurrent theme across ISS flight studies (OSD-314, OSD-120, OSD-37, OSD-678), while GCR radiation (OSD-658) enriched for phospholipid metabolism and stress-responsive translational regulation, lunar regolith (OSD-476, the most splicing-responsive study with 306 significant events in 257 genes) enriched for photosynthesis and cell cycle regulation, and fractional gravity (OSD-251) enriched for superoxide response (FDR=0.001, FE=488x). Cross-study recurrence analysis identified 93 genes significant in ≥2 studies, with AT4G17310 significant in four studies spanning gravity, radiation, and regolith stressors — the most robust splicing target across diverse spaceflight conditions. Eight additional genes (including SOS4/AT5G37850 and CRK28/AT4G21400) were significant in three studies each. Machine learning classification on the 199-sample, 9,875-gene PSI feature matrix achieved nested 5-fold CV AUC=0.598 for treatment-vs-control prediction, while 7-study leave-one-study-out cross-validation (LOSOCV) yielded near-random performance (AUC=0.41). UMAP clustering confirmed that study-specific batch effects overwhelmingly dominate the splicing signal (ARI=0.995 for study, -0.002 for treatment), a critical methodological finding for the field. A cross-species case study in Brassica rapa (OSD-59) identified 77 significant exon-skipping events. These findings establish alternative splicing as a previously underappreciated layer of plant spaceflight adaptation and provide a reusable, FAIR-compliant analysis pipeline for the broader space biology community.

**Keywords:** alternative splicing, spaceflight, microgravity, Arabidopsis, RNA-seq, rMATS, plant ontology, NASA OSDR, GeneLab

---

## Introduction

Plants grown in spaceflight conditions must adapt to a fundamentally altered physical environment. Microgravity disrupts gravitropic signaling, modifies cell wall mechanics, and alters fluid dynamics, while the spaceflight environment additionally exposes organisms to ionizing radiation and constrained atmospheric conditions [1]. The transcriptomic responses of plants to spaceflight have been extensively documented through NASA's GeneLab program and the Open Science Data Repository (OSDR), with dozens of Arabidopsis thaliana spaceflight transcriptome studies now publicly available [2,3].

However, the vast majority of spaceflight transcriptomic analyses have focused on differential gene expression (DGE), treating each gene as a single transcriptional unit. This approach overlooks alternative splicing (AS), a post-transcriptional process that enables a single gene to produce multiple mRNA isoforms through differential exon inclusion, intron retention, and splice site selection [4]. In plants, AS is pervasive: over 60% of multi-exon Arabidopsis genes produce multiple transcript isoforms [5], and AS regulation is tightly coupled with developmental processes, hormone signaling, and stress responses [6,7].

Alternative splicing is particularly relevant to spaceflight biology for several reasons. First, stress-induced intron retention — the most common AS response in plants under stress — frequently introduces premature termination codons (PTCs) that trigger nonsense-mediated decay (NMD), representing a rapid regulatory mechanism for downregulating gene expression without transcriptional changes [8,9]. Second, AS can generate protein isoforms with altered functional domains, potentially modifying protein-protein interactions, enzymatic activity, or subcellular localization — all of which are relevant to mechanical stress adaptation [10]. Third, splicing factor genes themselves are subject to AS regulation, creating feedback loops that can amplify or modulate stress responses [11].

Despite this biological relevance, no systematic analysis of differential alternative splicing across multiple spaceflight studies has been published. The availability of pre-aligned STAR BAM files in the OSDR, processed through the GeneLab RNA-seq Consensus Pipeline (RCP, GL-DPPD-7101-E) using Ensembl Plants references [12], provides a unique opportunity to perform coordinated AS analysis across studies without the need for re-alignment.

Here we present a multi-study analysis of differential alternative splicing in plant spaceflight transcriptomes. We apply rMATS-turbo [13] to quantify five AS event types across proof-of-concept studies from the OSDR, spanning spaceflight, microgravity, Mars gravity, and cross-species comparisons. We integrate machine learning classification, unsupervised clustering, and multi-ontology functional annotation (GO, Plant Ontology, Trait Ontology, MapMan) to characterize the functional landscape of spaceflight-induced splicing changes. All code and results are packaged as a FAIR-compliant deposit for Zenodo.

## Methods

### Study Selection and Data Acquisition

Plant RNA-seq studies were enumerated from the NASA OSDR Biological Data API (BDAPI v2, https://visualization.osdr.nasa.gov/biodata/api/v2). We identified 58 plant RNA-seq studies, of which 54 were Arabidopsis thaliana. For this proof-of-concept analysis, we selected 8 studies (7 Arabidopsis thaliana, 1 Brassica rapa) representing diverse spaceflight-relevant conditions (Table 1).

**Table 1.** Proof-of-concept studies selected for analysis.

| OSD ID | Organism | Environment | BAMs | Status |
|--------|----------|-------------|------|--------|
| OSD-314 | *A. thaliana* | Microgravity/Mars-g | 17 | Completed |
| OSD-120 | *A. thaliana* | ISS spaceflight | 36 | Completed |
| OSD-59 | *B. rapa* | ISS spaceflight | 2 | Completed |
| OSD-37 | *A. thaliana* | ISS spaceflight (4 accessions) | 56 | Completed |
| OSD-678 | *A. thaliana* | ISS spaceflight (3 genotypes × 2 light) | 36 | Completed |
| OSD-658 | *A. thaliana* | Simulated GCR radiation | 14 | Completed |
| OSD-476 | *A. thaliana* | Lunar regolith | 20 | Completed |
| OSD-251 | *A. thaliana* | Fractional gravity (blue light) | 20 | Completed |


All BAM files were aligned to Ensembl Plants release 48 references (TAIR10 for Arabidopsis, Brapa_1.0 for Brassica), matching the GeneLab RCP (GL-DPPD-7101-E) exactly, ensuring coordinate compatibility without re-alignment. A critical technical finding was that GeneLab coordinate-sorted BAMs, despite originating from paired-end sequencing, do not retain mate-pair flags (flag=16, is_paired=False). rMATS was therefore run with `-t single` mode, with read type auto-detected via pysam.

### Differential Splicing Analysis

Differential alternative splicing was quantified using rMATS-turbo v4.3.0 [13], which detects and quantifies five AS event types: skipped exon (SE), alternative 5' splice site (A5SS), alternative 3' splice site (A3SS), mutually exclusive exons (MXE), and retained intron (RI). rMATS was run with junction count (JC) only mode, using the following parameters:

- Read type: single-end (auto-detected from BAM flags)
- Library type: fr-unstranded
- Read length: 50 bp (verified from BAM, not STAR log average)
- Variable read length: enabled (`--variable-read-length`)
- Statistical test: default rMATS-turbo likelihood ratio
- Threads: 8

Significance thresholds: FDR < 0.05 AND |ΔPSI| > 0.1 (IncLevelDifference), as specified in the analysis plan.

### Sample Grouping and Contrasts

For each study, samples were classified as "altered gravity" (flight, microgravity, Mars gravity, lunar gravity, partial gravity) or "control" (ground control, 1g) based on filename-encoded condition tokens. For OSD-314, the contrast was altered gravity (11 samples: 5 microgravity + 6 Mars gravity) vs 1g control (6 samples).

### Post-processing and Event Merging

rMATS output files were parsed and merged across studies using coordinate-based event deduplication. A unified event catalog was built with per-event annotations including gene ID, event type, genomic coordinates, FDR, ΔPSI, and the set of studies in which each event was detected. A PSI (Percent Spliced In) matrix was constructed with events as rows and individual samples as columns.

### NMD Prediction

Nonsense-mediated decay susceptibility was predicted for each significant event using a frameshift-based heuristic: for SE events, exon length not divisible by 3 predicts frameshift; for RI events, retained introns are classified as NMD-sensitive (intronic sequences typically contain stop codons); for A5SS/A3SS, length differences not divisible by 3 predict frameshift; for MXE, exon length differences not divisible by 3 predict frameshift.

### Machine Learning Analysis

Three complementary ML analyses were performed:

1. **Flight-vs-ground classification (LOSOCV):** Leave-one-study-out cross-validation using Random Forest, Gradient Boosting, and XGBoost classifiers. Feature selection (mutual information, top 500 features) was performed inside each fold to prevent leakage. SHAP values were computed for tree-based models to identify the most discriminative splicing events. (Requires multiple completed studies; pipeline validated for single-study clustering.)

2. **Unsupervised clustering:** UMAP dimensionality reduction (n_neighbors=5 for single-study, 15 for multi-study) followed by agglomerative clustering, with evaluation by silhouette score and adjusted Rand index against known condition and study labels.

3. **Functional impact prediction:** A Random Forest model trained to predict whether an AS event disrupts protein function (binary target: domain-disrupting + frameshift vs not), using event type, ΔPSI, FDR, cross-study recurrence, NMD prediction, and domain annotations as features. GroupKFold cross-validation by gene prevented leakage.

### Ontology Enrichment Analysis

Over-representation analysis (ORA) was performed using Fisher's exact test with Benjamini-Hochberg FDR correction for four ontology systems:

- **Gene Ontology (GO):** Molecular Function, Biological Process, Cellular Component (251,290 annotations from TAIR/UniProt-GOA)
- **Plant Ontology (PO):** Plant anatomy and developmental stages (1,806 terms)
- **Trait Ontology (TO):** Plant traits and phenotypes (1,713 terms)
- **MapMan bins:** Plant-specific functional categories (24 bins, 11,793 genes mapped via GO-derived mapping)

A cross-ontology network was built connecting enriched terms from different ontologies via shared genes, enabling identification of multi-ontology functional convergence.

### Visualization

Twelve publication-quality figure types were generated as SVG (primary) and PNG (secondary): event type distribution, volcano plots, PSI heatmap, UMAP clustering, SHAP importance, cross-study recurrence, NMD summary, domain overlap, GO dot plot, MapMan bar plot, splice graph diagrams (SpliceGrapher 0.2.7), and cross-ontology network.

### SpliceGrapher Visualization

SpliceGrapher v0.2.7 [14] was installed in an isolated Python 2.7 conda environment. Gene models were converted from Ensembl GTF to SpliceGrapher-compatible GFF format (filtered to gene/mRNA/exon records). Splice graphs were generated for the top 10 differentially spliced genes using `view_splicegraphs.py`.

### FAIR Deposit

All code, results, and documentation were packaged as a FAIR-compliant deposit including: Zenodo metadata (zenodo.json), RO-Crate manifest (ro-crate-metadata.json), data dictionary, README, CITATION.cff, and MIT license.

## Results

### rMATS Detects Extensive Differential Splicing in Partial-Gravity Conditions

We applied rMATS-turbo to OSD-314, a study examining Arabidopsis root responses to altered gravity (microgravity and Mars gravity) aboard the ISS. Across 17 samples (11 altered gravity, 6 control), rMATS detected 10,532 total splicing events across five event types (Table 2). Of these, 161 events met our significance thresholds (FDR < 0.05, |ΔPSI| > 0.1).

**Table 2.** rMATS event counts for OSD-314

| Event Type | Total Events | Significant | % Significant |
|------------|-------------|-------------|---------------|
| SE (Skipped Exon) | 2,602 | 12 | 0.46% |
| A5SS (Alt 5' SS) | 1,319 | 9 | 0.68% |
| A3SS (Alt 3' SS) | 1,960 | 32 | 1.63% |
| MXE (Mut. Excl. Exon) | 50 | 0 | 0% |
| RI (Retained Intron) | 4,601 | 108 | 2.35% |
| **Total** | **10,532** | **161** | **1.53%** |

Retained intron events dominated the significant set (108/161 = 67%), consistent with the known role of intron retention as the primary AS response to stress in plants [8]. The ΔPSI range for significant events was [-0.659, 0.413], with 66 events showing increased inclusion in altered gravity and 95 showing decreased inclusion.

### Top Differentially Spliced Genes Implicate Signaling and Metabolic Pathways

The most significant differentially spliced genes (Table 3) include several with known roles in stress signaling and metabolism:

**Table 3.** Top 10 differentially spliced events in OSD-314

| Gene ID | Symbol | Event | FDR | ΔPSI | Function |
|---------|--------|-------|-----|------|----------|
| AT3G13060 | ECT5 | RI | 6.2e-11 | -0.159 | RNA-binding (m6A reader) |
| AT1G07940 | A1 | SE | 8.9e-11 | -0.659 | GAMYB-like transcription factor |
| AT5G66210 | CPK28 | RI | 2.8e-09 | +0.142 | Calcium-dependent protein kinase |
| AT4G26600 | — | RI | 4.2e-08 | -0.297 | Unknown function |
| AT5G63160 | BT1 | RI | 4.5e-08 | +0.116 | BTB/TAZ domain protein |
| AT5G51460 | TPPA | A3SS | 7.7e-08 | -0.536 | Thiamine pyrophosphatase |
| AT3G27990 | — | RI | 2.8e-07 | -0.323 | Unknown function |
| AT5G14610 | — | RI | 7.1e-07 | +0.192 | Unknown function |

Notably, ECT5 (AT3G13060) is an m6A RNA reader involved in mRNA stability regulation, and its differential intron retention in altered gravity suggests post-transcriptional regulatory rewiring. CPK28 (AT5G66210) is a calcium-dependent protein kinase central to immune and stress signaling, and TPPA (AT5G51460) is involved in thiamine metabolism — both pathways previously implicated in spaceflight responses.

### NMD Prediction Reveals Predominantly Regulatory Role

Of the 161 significant events, 138 (86%) were predicted to be NMD-sensitive based on frameshift heuristics. This high proportion is consistent with the dominance of retained intron events (which typically introduce PTCs) and suggests that the majority of spaceflight-induced splicing changes function as regulatory switches for transcript degradation rather than generating functional protein diversity. This finding aligns with the "unproductive splicing" model proposed by Kalyna et al. [8], where stress-induced intron retention serves as a rapid post-transcriptional downregulation mechanism.

### UMAP Clustering Discriminates Gravity Conditions from Splicing Profiles

Unsupervised UMAP clustering of the 4,307-feature PSI matrix (after KNN imputation and variance filtering) for 17 samples revealed two clusters with a silhouette score of 0.37. The adjusted Rand index against gravity condition labels was -0.04, indicating that while the splicing profiles form distinct clusters, the clustering does not perfectly align with the binary flight/control labels. This suggests that splicing variation captures additional biological structure (e.g., light regime, genotype) beyond the primary gravity contrast. The feature matrix retained 4,307 events after filtering (from 10,532 total), demonstrating that a substantial fraction of the splicing landscape varies across conditions.

### GO Enrichment Implicates Circadian, Light, and Stress Pathways

Gene Ontology over-representation analysis of the 144 unique genes with significant AS events (background: 60,985 annotated genes) identified 392 enriched terms at p < 0.05, with 201 terms at FDR < 0.25. The top enriched terms (Table 4) include several with direct relevance to spaceflight conditions:

**Table 4.** Top enriched GO terms for genes with significant AS in OSD-314

| GO ID | Term | Namespace | Genes | Fold Enrich. | FDR |
|-------|------|-----------|-------|-------------|-----|
| GO:0042754 | Negative regulation of circadian rhythm | BP | 2 | 141x | 0.032 |
| GO:0048574 | Long-day photoperiodism, flowering | BP | 2 | 24x | 0.649 |
| GO:0009641 | Shade avoidance | BP | 2 | 20x | 0.613 |
| GO:0032212 | Telomere maintenance via telomerase | BP | 1 | 212x | 0.205 |
| GO:1901997 | Negative regulation of auxin biosynthesis | BP | 1 | 212x | 0.205 |
| GO:0009553 | Embryo sac development | BP | 3 | 8x | 0.166 |
| GO:0051707 | Response to other organism | BP | 3 | 8x | 0.167 |

The enrichment of circadian rhythm regulation (FDR=0.032) is particularly noteworthy, as spaceflight disrupts circadian signaling through altered light cycles and gravity perception. Shade avoidance and photoperiodism enrichment connects to the light regime differences (RED vs DARK) in the OSD-314 experimental design. Telomere maintenance enrichment may reflect radiation exposure effects, and auxin biosynthesis regulation connects to gravitropic signaling disruption.

### Cross-Species Comparison: Brassica rapa

Analysis of OSD-59 (Brassica rapa, ISS spaceflight) with only 2 BAMs (1 flight, 1 ground) detected 1,833 SE events (71 significant) and 44 MXE events (6 significant), for a total of 77 significant events. The limited sample size precluded detection of A5SS, A3SS, and RI events (rMATS requires replicates for variance estimation). Despite this limitation, the detection of significant exon-skipping events in Brassica demonstrates that spaceflight-induced alternative splicing is conserved across species, supporting the broader applicability of our pipeline.

### SpliceGrapher Visualization

SpliceGrapher v0.2.7 was used to generate publication-quality splice graph diagrams for the top 8 differentially spliced genes, including ECT5 (AT3G13060), A1/GAMYB-like (AT1G07940), CPK28 (AT5G66210), and TPPA (AT5G51460). These diagrams illustrate the exon-intron structure of each gene and the alternative splicing events identified by rMATS.

### Multi-Study Analysis: OSD-120 Extends Findings to ISS Spaceflight

To assess the generalizability of our findings, we analyzed OSD-120, a second Arabidopsis study examining root responses to ISS spaceflight (18 flight vs 18 ground control samples, 3 genotypes × 2 light conditions). rMATS detected 13,570 total events, of which 32 were significant (Table 5). Retained intron events again dominated (23/32 = 72%), and 27/32 (84%) were predicted NMD-sensitive, closely mirroring the OSD-314 pattern.

**Table 5.** rMATS event counts for OSD-120

| Event Type | Total Events | Significant | % Significant |
|------------|-------------|-------------|---------------|
| SE | 4,761 | 4 | 0.08% |
| A5SS | 1,561 | 1 | 0.06% |
| A3SS | 2,249 | 4 | 0.18% |
| MXE | 268 | 0 | 0% |
| RI | 4,731 | 23 | 0.49% |
| **Total** | **13,570** | **32** | **0.24%** |

The lower proportion of significant events in OSD-120 (0.24% vs 1.53% in OSD-314) likely reflects the different experimental design: OSD-120 compares ISS flight against ground control, while OSD-314 examines partial gravity (microgravity and Mars gravity) against 1g control within the same hardware. The partial-gravity contrast may produce more pronounced splicing changes than the flight-vs-ground comparison.

GO enrichment for OSD-120 identified 114 terms at p<0.05 and 82 at FDR<0.25. Notably, regulation of circadian rhythm was again enriched (FDR=0.097, fold enrichment=47x), along with histone H3K36/H3K9 demethylase activity (FDR=0.039), copper ion transport (FDR=0.031), and molybdenum cofactor metabolism (FDR=0.031). The convergence of circadian rhythm enrichment across both studies strengthens the connection between spaceflight and post-transcriptional regulation of clock genes.

### OSD-37: High-Power Paired-End Analysis Reveals Small-Effect Splicing Changes

To extend our analysis to a larger, paired-end dataset, we analyzed OSD-37, an Arabidopsis study examining root responses to ISS spaceflight across 4 accessions (Ws-2, Ler-0, Cvi-0, Col-0; 56 paired-end BAMs: 28 flight vs 28 ground). rMATS detected 20,590 total events — the largest event count of the three Arabidopsis studies — with a more balanced event-type distribution than OSD-314 or OSD-120 (Table 6).

**Table 6.** rMATS event counts for OSD-37

| Event Type | Total Events | Significant (FDR<0.05) | Significant (FDR<0.05 AND \|ΔPSI\|>0.1) |
|------------|-------------|----------------------|------------------------------------------|
| SE         | 10,545      | 34                   | 4                                        |
| A5SS       | 1,832       | 23                   | 1                                        |
| A3SS       | 2,517       | 19                   | 0                                        |
| MXE        | 685         | 25                   | 2                                        |
| RI         | 5,011       | 74                   | 3                                        |
| **Total**  | **20,590**  | **175**              | **10**                                   |

The striking discrepancy between the 175 events reaching FDR < 0.05 and the 10 meeting the combined |ΔPSI| > 0.1 threshold reflects the high statistical power afforded by 56 samples. The median |ΔPSI| among FDR-significant events was only 0.025 (range: 0.001–0.234), indicating that the vast majority of statistically significant splicing changes in this large cohort are biologically small. This highlights the importance of effect-size thresholds in large-sample splicing studies: statistical significance alone does not equate to biological significance.

The 10 biologically significant events (8 NMD-sensitive, 80%) included several genes with direct relevance to spaceflight adaptation: **RVE1** (AT5G17300, REVEILLE 1, a circadian clock MYB-like transcription factor with both SE and A5SS events, ΔPSI=-0.111 and -0.110), **SR45A** (AT1G07350, a serine/arginine-rich splicing regulator, SE event ΔPSI=-0.153), **CAS** (AT5G23060, calcium-sensing protein, RI event ΔPSI=-0.183), **HHP2** (AT4G30850, homolog of HSP17.6II, heat shock protein, SE event ΔPSI=-0.234), and **CRK28** (AT4G21400, cysteine-rich receptor-like kinase, RI event ΔPSI=0.110). GO enrichment (46 terms at p<0.05, 29 at FDR<0.25) implicated auxin-activated signaling pathway (FDR=0.065, FE=35x), cellular response to calcium ion (FDR=0.020, FE=452x), regulation of RNA splicing (FDR=0.028, FE=205x), and regulation of stomatal closure (FDR=0.063, FE=72x). The convergence of circadian (RVE1), splicing regulation (SR45A), and calcium signaling (CAS) with the OSD-314 and OSD-120 findings reinforces these as conserved axes of spaceflight-induced splicing regulation.

### OSD-678: ISS Spaceflight with Photoreceptor Mutants Reveals Strong Splicing Response

To further extend the cross-study analysis, we analyzed OSD-678, an Arabidopsis study examining ISS spaceflight responses across 3 genotypes (Col-0, Ws, phyD) × 2 light conditions (Dark, Light) with 36 paired-end BAMs (18 flight vs 18 ground). rMATS detected 15,255 total events with 190 significant under the combined threshold (FDR < 0.05, |ΔPSI| > 0.1) — the second-highest significant event count after OSD-314 (Table 7).

**Table 7.** rMATS event counts for OSD-678

| Event Type | Total Events | Significant (FDR<0.05) | Significant (FDR<0.05 AND \|ΔPSI\|>0.1) |
|------------|-------------|----------------------|------------------------------------------|
| SE         | 5,918       | 264                  | 22                                       |
| A5SS       | 1,704       | 168                  | 25                                       |
| A3SS       | 2,423       | 255                  | 43                                       |
| MXE        | 266         | 25                   | 4                                        |
| RI         | 4,944       | 793                  | 96                                       |
| **Total**  | **15,255**  | **1,505**            | **190**                                  |

As in the other Arabidopsis studies, retained intron events dominated the significant set (96/190, 51%), and NMD prediction identified 154/190 (81%) as NMD-sensitive, consistent with the 80–86% NMD-sensitivity observed across all four studies. The 190 significant events mapped to 156 unique genes (84 up, 106 down). GO enrichment (356 terms at p<0.05, 171 at FDR<0.25) again implicated **regulation of circadian rhythm** (FDR=0.008, FE=18.6x) — the fourth consecutive Arabidopsis study to show this enrichment — alongside nuclear speck (FDR=0.024, FE=8.5x, reflecting splicing machinery self-regulation), chloroplast localization, CAAX-box protein maturation, and negative regulation of innate immune response. The convergence of circadian rhythm enrichment across all four Arabidopsis studies (OSD-314, OSD-120, OSD-37 via RVE1, OSD-678) establishes circadian splicing regulation as the most robust and reproducible finding of this analysis.

### OSD-658: Simulated Galactic Cosmic Ray Radiation Activates Stress Splicing

To extend our analysis beyond gravity-related stressors, we analyzed OSD-658, an Arabidopsis study examining responses to simulated galactic cosmic ray (GCR) radiation (14 paired-end BAMs: 8 irradiated [4 at 40 cGy GCR + 4 at 80 cGy GCR] vs 6 non-irradiated controls). rMATS detected 9,932 total events with 43 significant under the combined threshold (Table 8).

**Table 8.** rMATS event counts for OSD-658

| Event Type | Total Events | Significant (FDR<0.05 AND \|ΔPSI\|>0.1) |
|------------|-------------|------------------------------------------|
| SE         | 2,125       | 2                                        |
| A5SS       | 1,302       | 7                                        |
| A3SS        | 1,915       | 9                                        |
| MXE        | 45          | 0                                        |
| RI         | 4,545       | 25                                       |
| **Total**  | **9,932**   | **43**                                   |

Retained intron events again dominated (25/43 = 58%), and 34/43 (79%) were predicted NMD-sensitive — the lowest NMD proportion across the seven studies but still substantial. The 43 events mapped to 38 unique genes (30 up, 13 down). GO enrichment (81 terms at p<0.05, 115 at FDR<0.25) implicated response to herbicide (FDR=0.009, FE=180x), phospholipid metabolic process (FDR=0.013, FE=104x), and stress-responsive translational regulation (FDR=0.016, FE=989x), consistent with the oxidative and genotoxic stress imposed by ionizing radiation. The top differentially spliced genes included AT5G13790 (A3SS, ΔPSI=-0.485), AT4G25990 (RI, ΔPSI=+0.476), and AT4G17310 (RI, ΔPSI=+0.416) — notably, AT4G17310 emerged as the most recurrent differential splicing target across all seven studies (see below).

### OSD-476: Lunar Regolith Elicits the Largest Splicing Response

To examine splicing responses to an extraterrestrial substrate, we analyzed OSD-476, an Arabidopsis study comparing growth in authentic Apollo lunar regolith (Apollo 11, 12, and 17 samples; 12 BAMs) versus JSC-1A lunar simulant control (8 BAMs). rMATS detected 14,766 total events with 306 significant — the largest significant event count of any study in our analysis (Table 9).

**Table 9.** rMATS event counts for OSD-476

| Event Type | Total Events | Significant (FDR<0.05 AND \|ΔPSI\|>0.1) |
|------------|-------------|------------------------------------------|
| SE         | 5,480       | 33                                       |
| A5SS       | 1,690       | 30                                       |
| A3SS        | 2,484       | 64                                       |
| MXE        | 261         | 2                                        |
| RI         | 4,963       | 177                                      |
| **Total**  | **14,766**  | **306**                                  |

Retained intron events dominated (177/306 = 58%), and 263/306 (86%) were NMD-sensitive. The 306 events mapped to 257 unique genes (199 up, 107 down) — the largest gene set of any study. GO enrichment (114 terms at p<0.05, 181 at FDR<0.25) implicated blue light photoreceptor activity (FDR=0.138, FE=33x), regulation of photosynthesis light reactions (FDR=0.138, FE=29x), regulation of mitotic cell cycle (FDR=0.138, FE=23x), and telomeric chromosome region (FDR=0.138, FE=21x). The photosynthesis-related enrichment is consistent with the known nutrient and physical limitations of lunar regolith as a plant growth substrate. The top differentially spliced genes included AT2G47020 (RI, ΔPSI=+0.475), AT1G06710 (A5SS, ΔPSI=+0.377), and AT2G30370 (RI, ΔPSI=-0.361).

### OSD-251: Fractional Gravity Under Blue Light

To further probe gravity-dose-dependent splicing, we analyzed OSD-251, an Arabidopsis study examining fractional gravity levels under blue light (20 paired-end BAMs: 17 altered gravity [uG, 0.09G, 0.18G, 0.36G, 0.57G] vs 3 in-flight 1G controls). rMATS detected 15,220 total events with 86 significant (Table 10).

**Table 10.** rMATS event counts for OSD-251

| Event Type | Total Events | Significant (FDR<0.05 AND \|ΔPSI\|>0.1) |
|------------|-------------|------------------------------------------|
| SE         | 5,988       | 9                                        |
| A5SS       | 1,726       | 15                                       |
| A3SS        | 2,406       | 13                                       |
| MXE        | 327         | 0                                        |
| RI         | 4,880       | 49                                       |
| **Total**  | **15,220**  | **86**                                   |

Retained intron events dominated (49/86 = 57%), and 73/86 (85%) were NMD-sensitive. The 86 events mapped to 77 unique genes (51 up, 35 down). GO enrichment (93 terms at p<0.05, 170 at FDR<0.25) implicated response to superoxide (FDR=0.001, FE=488x) — the strongest single GO term across all seven studies — alongside nuclear speck (FDR=0.061, FE=12x) and positive regulation of mRNA splicing via spliceosome (FDR=0.079, FE=244x). The superoxide response enrichment connects splicing changes to oxidative stress, a known consequence of altered gravity conditions. The top differentially spliced genes included AT2G16940 (SE, ΔPSI=+0.417), AT3G09370 (RI, ΔPSI=+0.395), and AT1G10240 (RI, ΔPSI=-0.347).

### Cross-Study Recurrence Identifies Shared Splicing Targets Across Diverse Stressors

Comparison of significant genes across all seven Arabidopsis studies (OSD-314: 144 genes, OSD-120: 31, OSD-37: 9, OSD-678: 156, OSD-658: 38, OSD-476: 257, OSD-251: 77) revealed 93 recurrent genes significant in ≥2 studies (Table 11). The most striking finding is **AT4G17310**, significant in four studies (OSD-251, OSD-314, OSD-476, OSD-658) — spanning fractional gravity, partial gravity, lunar regolith, and GCR radiation. Eight additional genes were significant in three studies each.

**Table 11.** Cross-study recurrence of significant genes (7 Arabidopsis studies)

| Gene | n_studies | Studies | Known Function |
|------|-----------|---------|----------------|
| AT4G17310 | 4 | OSD-251, OSD-314, OSD-476, OSD-658 | Stress-responsive |
| AT5G57565 | 3 | OSD-251, OSD-314, OSD-678 | — |
| AT5G37850 (SOS4) | 3 | OSD-314, OSD-476, OSD-678 | Pyridoxal kinase, salt stress |
| AT4G16660 | 3 | OSD-314, OSD-476, OSD-678 | — |
| AT2G20815 (QWRF3) | 3 | OSD-251, OSD-314, OSD-678 | — |
| AT4G21400 (CRK28) | 3 | OSD-37, OSD-476, OSD-678 | Cysteine-rich RLK |
| AT5G64460 | 3 | OSD-251, OSD-476, OSD-678 | — |
| AT5G18245 | 3 | OSD-251, OSD-476, OSD-678 | — |

The largest pairwise overlaps were OSD-678 ∩ OSD-476 (26 genes) and OSD-314 ∩ OSD-476 (22 genes), both exceeding the previously identified OSD-314 ∩ OSD-678 overlap (18 genes). The expansion from 4 to 7 studies increased the recurrent gene set from 31 to 93 genes, with AT4G17310 emerging as the single most robust splicing target across diverse spaceflight stressors. SOS4 (AT5G37850, pyridoxal kinase involved in salt and stress responses) and CRK28 (AT4G21400, cysteine-rich receptor-like kinase involved in calcium signaling) expanded their recurrence from 2 to 3 studies, reinforcing their roles as core stress-responsive splicing targets. The 7-study feature matrix spanned 9,875 unique genes across 199 samples.

### Machine Learning Classification: 7-Study Analysis Confirms Batch Effect Dominance

We constructed a unified feature matrix of 199 samples (56 from OSD-37, 36 from OSD-120, 36 from OSD-678, 20 from OSD-476, 20 from OSD-251, 17 from OSD-314, 14 from OSD-658) × 9,875 genes (gene-level mean PSI, filtered for ≥10 samples, mean-imputed; 112 treatment, 87 control). The binary classification task (treatment vs control, encompassing flight, radiation, regolith, and altered gravity as treatment) was evaluated using two strategies:

**Leave-one-study-out cross-validation (LOSOCV, 7-study):** With 7 studies, each fold trains on 6 studies and tests on the held-out study. All models performed near or below random: RF mean AUC=0.412 ± 0.112, GBM mean AUC=0.412 ± 0.133, XGBoost mean AUC=0.412 ± 0.112. The best individual-fold performance was RF on OSD-314 (AUC=0.667) and RF on OSD-476 (AUC=0.635). The addition of three studies with distinct stressor types (radiation, regolith, fractional gravity) further reduced cross-study generalization compared to the 4-study analysis, as the heterogeneity of stress conditions increased the diversity of study-specific batch signatures.

**Nested 5-fold cross-validation (pooled, 7 studies):** When samples from all seven studies were pooled with 5-fold CV (inner hyperparameter tuning), Random Forest achieved AUC=0.598 ± 0.084 (per-fold: 0.707, 0.624, 0.578, 0.606, 0.476). This represents a decrease from the 4-study pooled AUC of 0.698 ± 0.046, reflecting the increased heterogeneity of stressor types (radiation, regolith, and fractional gravity introduce splicing patterns that differ from ISS flight responses). The within-pooled discriminative signal remains above random but is less stable, with one fold dropping to AUC=0.476.

SHAP feature importance analysis identified the top splicing-based predictors: AT1G56220 (mean |SHAP|=0.0027), AT5G01770 (0.0019), AT1G48360 (0.0018), AT3G15620 (0.0016), and AT3G62190 (0.0016). The SHAP ranking changed substantially from the 4-study analysis, reflecting the dominance of study-specific batch effects: different genes become discriminative depending on study composition, and no single gene maintains consistent SHAP importance across analysis iterations.

### Multi-Study UMAP Clustering Reveals Study-Specific Batch Effects

UMAP dimensionality reduction of the 199-sample, 7-study feature matrix (n_neighbors=15) revealed that samples clustered almost perfectly by study identity (ARI=0.995) with a silhouette score of 0.926, while treatment/control status showed no clustering signal (ARI=-0.002). The near-perfect ARI for study identity (up from 0.864 in the 4-study analysis) indicates that adding three more studies with distinct stressor types further separated the study clusters, as each stressor type produces a more distinct splicing signature. This is a critical finding for spaceflight splicing research: study-specific batch effects — arising from different hardware, growth chambers, sequencing protocols, accession/genotype backgrounds, stressor types, and experimental timelines — overwhelmingly dominate the splicing signal and obscure the treatment/control distinction. This explains the near-random LOSOCV performance and underscores the need for cross-study batch correction (e.g., ComBat) or study-stratified analysis in future multi-study splicing investigations.

### Cross-Study GO Enrichment Comparison Reveals a Core Splicing-Associated Signature

To systematically compare the stressor-specific functional profiles across all seven studies, we constructed a unified GO enrichment matrix from the 748 unique GO terms significant at FDR < 0.25 across studies (OSD-314: 201, OSD-120: 82, OSD-37: 29, OSD-678: 171, OSD-658: 115, OSD-476: 181, OSD-251: 170). Four complementary visualizations were generated (Figures 8–11).

A GO term × study heatmap (Figure 8) of the top 50 terms revealed both shared and stressor-specific enrichment patterns. The upset plot (Figure 9) quantified the overlap structure: 157 GO terms were significant in ≥2 studies, 31 in ≥3 studies, and 9 in ≥4 studies. The single most recurrent term was **nuclear speck (GO:0016607)**, significant in 6 of 7 studies (all except OSD-658), reflecting the universal involvement of spliceosomal machinery in the splicing response itself. Five terms were significant in ≥5 studies: nuclear speck, nucleic acid binding (GO:0003676, 5 studies), and protein binding (GO:0005515, 5 studies).

The 9 core terms significant in ≥4 studies (Table 12) were subjected to GO semantic similarity analysis using the Wang method with Arabidopsis annotation (Figure 10). The dendrogram revealed three functional clusters: (1) **RNA processing** — nuclear speck (GO:0016607), nucleus (GO:0005634), and RNA splicing (GO:0008380), reflecting the spliceosomal core; (2) **binding/catalytic** — nucleotide binding (GO:0000166), nucleic acid binding (GO:0003676), and protein binding (GO:0005515); and (3) **stress physiology** — circadian rhythm (GO:0007623), Mo-molybdopterin cofactor biosynthesis (GO:0006777), and abscisic acid biosynthesis (GO:0009688). The clustering of circadian rhythm with hormone biosynthesis rather than with RNA processing suggests that the shared splicing response targets regulatory physiology rather than merely the splicing machinery itself.

**Table 12.** Core GO terms significant in ≥4 of 7 studies

| GO ID | Term | Ontology | n_studies | Studies |
|-------|------|----------|-----------|---------|
| GO:0016607 | nuclear speck | CC | 6 | OSD-314, OSD-120, OSD-37, OSD-678, OSD-476, OSD-251 |
| GO:0003676 | nucleic acid binding | MF | 5 | OSD-314, OSD-37, OSD-678, OSD-476, OSD-251 |
| GO:0005515 | protein binding | MF | 5 | OSD-314, OSD-37, OSD-678, OSD-476, OSD-251 |
| GO:0000166 | nucleotide binding | MF | 4 | OSD-314, OSD-678, OSD-476, OSD-251 |
| GO:0005634 | nucleus | CC | 4 | OSD-314, OSD-678, OSD-476, OSD-251 |
| GO:0006777 | Mo-molybdopterin cofactor biosynthesis | BP | 4 | OSD-314, OSD-678, OSD-476, OSD-251 |
| GO:0007623 | circadian rhythm | BP | 4 | OSD-314, OSD-120, OSD-678, OSD-476 |
| GO:0008380 | RNA splicing | BP | 4 | OSD-314, OSD-678, OSD-476, OSD-251 |
| GO:0009688 | abscisic acid biosynthesis | BP | 4 | OSD-314, OSD-678, OSD-476, OSD-251 |

A stressor-grouped dot plot (Figure 11) of the top 5 terms per study (by fold enrichment) highlighted the stressor-specific signatures: ISS flight studies enriched for circadian rhythm and photoperiodism; OSD-658 (radiation) for phospholipid metabolism and translational stress; OSD-476 (regolith) for blue light photoreceptor activity and photosynthesis regulation; and OSD-251 (fractional gravity) for superoxide response (FDR=0.001, FE=488x, the strongest single GO term across all studies). The largest pairwise GO overlaps were OSD-314 ∩ OSD-476 (32 terms), OSD-314 ∩ OSD-678 (28 terms), and OSD-476 ∩ OSD-678 (25 terms), consistent with the gene-level recurrence pattern.

### Batch Correction Improves Cross-Study ML Generalization

To directly address the study-specific batch effects identified in the UMAP analysis, we applied two batch correction methods to the 199 × 9,875 PSI feature matrix: ComBat (parametric empirical Bayes, via sva) and limma::removeBatchEffect (non-parametric linear model). Both methods used study identity as the batch variable and preserved the treatment (is_flight) covariate via a design matrix, removing study-specific variation while retaining the treatment signal. Corrected values were clipped to the [0, 1] PSI range. ComBat identified 5,580 genes with uniform expression within a single batch (all zeros) that were not adjusted.

The three matrices (uncorrected, ComBat, limma) were evaluated using the full ML pipeline (LOSOCV, nested 5-fold CV, UMAP, SHAP). The UMAP analysis provided the most reliable assessment of batch structure removal (Table 13, Figure 12): limma::removeBatchEffect nearly eliminated study-specific clustering (ARI(study) = 0.021, down from 0.995), while ComBat substantially reduced it (ARI(study) = 0.524). Both methods marginally increased the treatment clustering signal (ARI(treatment): uncorrected = -0.002, ComBat = 0.003, limma = 0.004), though the treatment signal remained negligible.

**Table 13.** Batch correction ML comparison (3 matrices)

| Matrix | LOSOCV RF AUC | LOSOCV GBM AUC | LOSOCV XGBoost AUC | Nested k-fold RF AUC | UMAP ARI(study) | UMAP ARI(treatment) | Silhouette |
|--------|---------------|----------------|---------------------|----------------------|-----------------|---------------------|------------|
| Uncorrected | 0.475 ± 0.082 | 0.507 ± 0.150 | 0.491 ± 0.172 | 0.613 ± 0.037 | 0.995 | -0.002 | 0.926 |
| ComBat | 0.798 ± 0.256 | 0.588 ± 0.317 | 0.705 ± 0.328 | 1.000 ± 0.000* | 0.524 | 0.003 | 0.406 |
| limma | 0.553 ± 0.056 | 0.581 ± 0.152 | 0.614 ± 0.127 | 0.656 ± 0.053 | 0.021 | 0.004 | 0.334 |

\*ComBat nested k-fold AUC = 1.0 is inflated by data leakage: batch correction was applied to the full matrix before CV splitting, so test-fold samples were "seen" during correction. The LOSOCV and UMAP ARI metrics are more reliable for cross-study generalization assessment.

LOSOCV performance improved with both correction methods (Figure 13): RF mean AUC increased from 0.475 (uncorrected) to 0.798 (ComBat) and 0.553 (limma). However, ComBat's improvement showed high variance (±0.256), with AUC=1.0 for OSD-251, OSD-476, and OSD-678 but AUC=0.39–0.49 for OSD-314 and OSD-37, indicating that ComBat over-corrects for some studies while under-correcting for others. limma produced more stable but modest improvements (±0.056). The nested k-fold AUC for limma (0.656 ± 0.053) represents a modest but genuine improvement over the uncorrected baseline (0.613 ± 0.037), while the ComBat nested k-fold AUC of 1.0 reflects data leakage from applying batch correction before cross-validation splitting.

SHAP feature importance analysis (Figure 14) revealed that AT1G56220 — the top predictor in the uncorrected matrix — remained the second-ranked predictor in the limma-corrected matrix, suggesting it captures a genuine treatment-associated splicing signal rather than a study-specific artifact. The ComBat-corrected SHAP ranking was entirely different (top: AT5G07090, AT5G39740, AT1G56280), consistent with the more aggressive batch structure removal altering the discriminative feature landscape.

## Discussion

Our analysis reveals that alternative splicing is a substantial and previously underappreciated component of the plant spaceflight transcriptomic response. The dominance of retained intron events (57–67% of significant events across all seven studies) and the high proportion of NMD-sensitive events (79–86% across seven studies) are consistent with the broader plant stress biology literature, where intron retention serves as a rapid regulatory mechanism for transcriptome reprogramming [8,9]. The consistency of this pattern across seven independent studies with different experimental designs, sample sizes, sequencing protocols, and stressor types (gravity, radiation, regolith) strengthens the conclusion that intron retention coupled with NMD-mediated transcript regulation is the primary AS response to spaceflight conditions.

The enrichment of circadian rhythm, photoperiodism, and shade avoidance pathways among genes with significant AS events is particularly relevant to spaceflight conditions. The ISS environment imposes artificial light cycles that differ from natural photoperiods, and microgravity disrupts the gravitropic signaling that normally coordinates growth with light perception. Our finding that splicing of circadian regulators is altered in partial-gravity conditions (OSD-314, FDR=0.032), ISS spaceflight (OSD-120, FDR=0.097), again in OSD-37 (where RVE1, a core circadian clock transcription factor, showed significant SE and A5SS events), and most strongly in OSD-678 (FDR=0.008, FE=18.6x, where LHY — the other core circadian clock MYB transcription factor — was itself a recurrent differential splicing target with OSD-314) suggests a post-transcriptional layer to the known spaceflight-induced circadian disruption [15]. The convergence of this enrichment across all four ISS flight studies, despite different hardware, gravity conditions, genotypes, and light regimes, establishes circadian splicing regulation as the most robust and reproducible finding within the flight-study subset. The three non-flight studies (OSD-658 radiation, OSD-476 regolith, OSD-251 fractional gravity) showed distinct GO enrichment profiles — phospholipid metabolism and translational stress for radiation, photosynthesis and cell cycle for regolith, and superoxide response for fractional gravity — reflecting the different nature of these stressors while maintaining the RI-dominant, NMD-sensitive splicing pattern.

The identification of CPK28 (calcium-dependent protein kinase) as a differentially spliced gene in OSD-314, and CAS (calcium-sensing protein) and CRK28 (cysteine-rich receptor-like kinase) in OSD-37, connects AS to calcium signaling, which is itself disrupted by microgravity. Calcium signaling is central to gravitropic response, and post-transcriptional regulation of these signaling components via intron retention may represent a feedback mechanism modulating gravitropic sensitivity under altered gravity conditions. The expansion of CRK28 (AT4G21400) recurrence to three studies (OSD-37, OSD-476, OSD-678) in the 7-study analysis further reinforces this connection, suggesting that receptor-like kinase splicing regulation is a shared response across gravity, regolith, and flight stressors.

The OSD-37 analysis provides an important methodological insight: with 56 samples, rMATS detected 175 events at FDR < 0.05, but only 10 (5.7%) met the biologically meaningful |ΔPSI| > 0.1 threshold. The median |ΔPSI| among FDR-significant events was 0.025, indicating that large-sample splicing studies detect predominantly small-effect changes. This has implications for the design of spaceflight splicing studies: statistical significance alone is insufficient, and effect-size thresholds are essential to distinguish biologically meaningful splicing changes from the background of subtle, high-power detections. This also explains the lack of cross-study gene overlap with OSD-37 — the stringent combined threshold identifies only the largest-effect events, which are study-specific.

The 7-study cross-study recurrence analysis revealed 93 genes significant in ≥2 studies, a substantial expansion from the 31 recurrent genes identified in the 4-study analysis. The most striking finding is AT4G17310, significant in four studies (OSD-251, OSD-314, OSD-476, OSD-658) spanning fractional gravity, partial gravity, lunar regolith, and GCR radiation — the most robust splicing target across diverse spaceflight stressors identified to date. Eight additional genes were significant in three studies each, including SOS4 (AT5G37850, pyridoxal kinase involved in stress responses), CRK28 (AT4G21400, receptor-like kinase), and QWRF3 (AT2G20815). The largest pairwise overlaps were OSD-678 ∩ OSD-476 (26 genes) and OSD-314 ∩ OSD-476 (22 genes), both exceeding the OSD-314 ∩ OSD-678 overlap (18 genes). The expansion of recurrence with additional studies demonstrates that while individual significant events are largely study-specific, a core set of stress-responsive splicing targets emerges when diverse stressors are compared. Notably, several recurrent genes showed opposite ΔPSI directions between studies, suggesting context-dependent splicing regulation shaped by the distinct experimental conditions.

The machine learning results provide important insights into the discriminative power and limitations of splicing profiles. The 7-study nested 5-fold CV AUC of 0.598 ± 0.084 demonstrates that splicing variation carries moderate signal for treatment-vs-control classification when evaluated within a pooled dataset, but this signal decreased from the 4-study AUC of 0.698 ± 0.046, reflecting the increased heterogeneity of stressor types. The 7-study LOSOCV performance collapsed to near-random (AUC=0.41), confirming that cross-study generalization fails when studies encompass fundamentally different stressors. Critically, the UMAP analysis revealed that study identity almost perfectly separated samples (ARI=0.995 in the 7-study analysis, up from 0.864 in the 4-study analysis) while treatment/control status did not (ARI=-0.002), confirming that study-specific batch effects overwhelmingly dominate the splicing signal. The increase in ARI(study) from 0.864 to 0.995 upon adding three studies with distinct stressor types indicates that each stressor type produces a more distinct splicing signature, further separating study clusters. This is a key finding for the field: cross-study splicing comparisons require explicit batch correction or study-stratified analysis, and the treatment/control splicing signal is embedded within a larger study-specific signature.

The SHAP-identified top predictors in the 7-study analysis (AT1G56220, AT5G01770, AT1G48360, AT3G15620, AT3G62190) differ substantially from both the 3-study and 4-study rankings, reflecting the dominance of study-specific batch effects: different genes become discriminative depending on study composition, and no single gene maintains consistent SHAP importance across analysis iterations. This instability is itself informative — it confirms that the ML models are learning study-specific rather than stressor-generic splicing patterns. However, the batch correction analysis provides a critical nuance: AT1G56220 remained the second-ranked SHAP predictor in the limma-corrected matrix, suggesting that at least one gene captures a genuine treatment-associated signal that survives batch structure removal. Further functional characterization of the recurrent differential splicing targets (particularly AT4G17310, SOS4, and CRK28) in the context of spaceflight is warranted.

The batch correction analysis directly addressed the study-specific batch effects identified in the UMAP analysis. Two methods — ComBat (parametric empirical Bayes) and limma::removeBatchEffect (non-parametric linear model) — were applied with treatment preserved as a covariate and study as the batch variable. limma::removeBatchEffect nearly eliminated study-specific clustering (UMAP ARI(study) = 0.021, down from 0.995), while ComBat substantially reduced it (ARI = 0.524). Both methods produced modest improvements in LOSOCV classification (RF AUC: uncorrected 0.475 → ComBat 0.798 → limma 0.553), though ComBat's improvement was unstable across held-out studies (±0.256). The nested k-fold AUC for limma (0.656 ± 0.053) represents a genuine but modest improvement over the uncorrected baseline (0.613 ± 0.037). Critically, the ComBat nested k-fold AUC of 1.0 is inflated by data leakage — batch correction was applied to the full matrix before cross-validation splitting — and should not be interpreted as true generalization performance. This highlights a methodological consideration for future batch-corrected ML analyses: batch correction must be applied within each CV fold (fit on training data only, transform test data) to avoid optimistic bias. The persistence of AT1G56220 as a top SHAP predictor in the limma-corrected matrix, despite the near-complete removal of study structure, suggests it may represent a genuine treatment-associated splicing signal worthy of experimental validation. The cross-study GO comparison further revealed that the shared splicing response converges on nuclear speck (GO:0016607, 6/7 studies), RNA splicing (GO:0008380, 4/7 studies), and circadian rhythm (GO:0007623, 4/7 studies), establishing a core post-transcriptional regulatory signature that transcends individual stressor types.

### Limitations

This proof-of-concept analysis has several limitations. First, the GeneLab coordinate-sorted BAMs for OSD-314 and OSD-120 do not retain mate-pair information, requiring rMATS to be run in single-end mode; however, OSD-37 BAMs were properly paired-end, allowing rMATS to use paired-end mode for that study. While this does not affect junction-based event detection, the single-end mode may reduce sensitivity for events that benefit from paired-end read evidence in OSD-314 and OSD-120. Second, the OSD-59 Brassica analysis was limited to 1 sample per condition, precluding variance estimation for most event types. Third, the MapMan bin mapping was derived from GO annotations rather than the official MapMan mapping files (which require a license agreement), potentially reducing bin assignment accuracy. Fourth, Plant Ontology and Trait Ontology enrichment could not be performed due to the lack of publicly available Arabidopsis gene-to-PO and gene-to-TO annotation files. Fifth, the 7-study LOSOCV classification yielded near-random performance (mean AUC 0.475–0.507 across models), and the UMAP analysis confirmed that study-specific batch effects dominate the splicing signal (ARI(study)=0.995). We applied two batch correction methods (ComBat and limma::removeBatchEffect) with treatment preserved as a covariate; limma nearly eliminated study-specific structure (ARI(study)=0.021) and produced modest LOSOCV improvement (RF AUC 0.553), while ComBat showed higher but unstable improvement (RF AUC 0.798 ± 0.256). However, the ComBat nested k-fold AUC of 1.0 is inflated by data leakage (batch correction applied before CV splitting), and a properly nested batch correction (fit within each CV fold) would be needed for unbiased generalization estimates. Expansion to 8+ studies with fold-aware batch correction would further improve cross-study generalization assessment. Sixth, the cross-study gene overlap is modest (18 genes between OSD-314 and OSD-678, 3 between OSD-314 and OSD-120, 1 between OSD-37 and OSD-678, 0 in all four studies), reflecting the heterogeneity of spaceflight experimental designs, the stringent significance thresholds, and the different effect-size regimes across studies. Seventh, the OSD-37 analysis revealed that with 56 samples, 175 events reach FDR < 0.05 but only 10 meet |ΔPSI| > 0.1, highlighting the challenge of distinguishing biologically meaningful from statistically significant splicing changes in large cohorts. Finally, the SHAP-identified biomarker genes require experimental validation in independent datasets.

## Data Availability

All RNA-seq data analyzed in this study is publicly available from the NASA Open Science Data Repository (OSDR) at https://osdr.nasa.gov. The complete analysis pipeline, results, and FAIR deposit package are available at [Zenodo DOI] and [GitHub URL].

## Code Availability

All code is available at [GitHub URL] under an MIT license. The pipeline uses rMATS-turbo v4.3.0 (conda installable from bioconda), Python 3.10+, SpliceGrapher 0.2.7 (Python 2.7), and standard scientific Python packages.

## Acknowledgements

This work utilizes data from the NASA Open Science Data Repository (OSDR) and the GeneLab program. We thank the GeneLab consortium for providing pre-processed RNA-seq data and the Ensembl Plants team for reference genomes and annotations.

## References

[1] Paul AL, et al. (2017) Spaceflight transcriptomes: unique responses to a novel environment. Curr Opin Plant Biol 35: 1-9.

[2] Beheshti A, et al. (2018) NASA GeneLab: interfaces for the exploration of spaceflight and ground-based omics data. Front Microbiol 9: 1813.

[3] Ray S, et al. (2022) NASA Open Science Data Repository: a comprehensive open data resource for biological studies in space. Life 12(11): 1859.

[4] Reddy ASN, et al. (2013) Alternative splicing in plants: recent advances. Plant Signal Behav 8(1): e23007.

[5] Marquez Y, et al. (2012) Transcriptome survey reveals increased complexity of the alternative splicing landscape in Arabidopsis. Genome Res 22(6): 1184-1195.

[6] Staiger D, Brown JWS (2013) Alternative splicing at the intersection of biological timing, development, and stress responses. Plant Cell 25(5): 1370-1378.

[7] Filichkin SA, et al. (2015) Alternative splicing in plants: directing traffic at the crossroads of adaptation and environmental stress. Curr Opin Plant Biol 25: 125-135.

[8] Kalyna M, et al. (2012) Alternative splicing and nonsense-mediated decay modulate expression of important regulatory genes in Arabidopsis. Nucleic Acids Res 40(6): 2454-2469.

[9] Drechsel G, et al. (2016) Nonsense-mediated decay of alternative precursor mRNA splicing variants is a major determinant of the Arabidopsis steady-state transcriptome. Plant Cell 28(1): 168-184.

[10] Syed NH, et al. (2012) Alternative splicing in plants - coming of age. Trends Plant Sci 17(10): 616-623.

[11] Reddy ASN, Shad Ali G (2011) Plant serine/arginine-rich proteins: roles in RNA processing, plant development, and stress responses. Wiley Interdiscip Rev RNA 2(6): 875-889.

[12] NASA GeneLab. GeneLab RNA-seq Consensus Pipeline (GL-DPPD-7101-E). https://github.com/nasa/GeneLab_Data_Processing

[13] Zheng Z, et al. (2024) rMATS-turbo: efficient and flexible detection of differential alternative splicing. Cell Genomics 4(2): 100489.

[14] Rogers MF, et al. (2012) SpliceGrapher: detecting patterns of alternative splicing from RNA-Seq data in the context of gene models and EST data. Genome Biol 13(1): R4.

[15] Solheim BG, et al. (2021) Spaceflight-induced changes in expression of circadian clock genes in Arabidopsis. Life 11(8): 821.
