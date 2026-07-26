# Plan: Plant Alternative Splicing Analysis System for NASA OSDR

## System to pull differential alternative splicing from plant RNA-seq studies in NASA OSDR, model splice variants with SpliceGrapher, apply ML + plant ontology analysis, and produce a FAIR Zenodo deposit + npj Microgravity manuscript.

---

## 1. Summary

Build a complete, reproducible computational system that retrieves plant RNA-seq data from the NASA Open Science Data Repository (OSDR), quantifies differential alternative splicing (AS) across spaceflight and space-analog environments, visualizes splice models with SpliceGrapher, applies multi-pronged machine learning and plant ontology enrichment to connect transcript-level changes to molecular/cellular/physiological adaptation, and packages everything as a FAIR Zenodo deposit with a submission-ready npj Microgravity manuscript (Markdown + LaTeX).

**Proof-of-concept execution**: The pipeline is designed for all 54 Arabidopsis RNA-seq studies with pre-aligned BAMs in OSDR, but will be run on a representative 6-8 study subset spanning spaceflight, microgravity, radiation, low pressure, and light treatments to generate real results for the manuscript. The full 54-study run is executable via the same code on a cluster.

---

## 2. Confirmed Decisions

| Decision | Choice |
|---|---|
| Study scope | All 54 Arabidopsis RNA-seq studies with BAMs (primary); 3 crop species (Brachypodium OSD-375, Brassica OSD-59, tomato OSD-767) as comparative case studies. Wheat OSD-622 excluded (microarray, not RNA-seq). |
| SpliceGrapher | Run in isolated Python 2.7 conda environment for splice graph visualization; rMATS provides quantitative differential splicing calls. |
| AS quantification tool | rMATS (event-based PSI, 5 event types: SE, A5SS, A3SS, MXE, RI, robust statistics). |
| ML target | All three: (1) flight-vs-ground classification, (2) unsupervised clustering, (3) functional impact prediction. |
| Reference annotation | TAIR10 assembly + Araport11 annotation via Ensembl Plants release 48 — exactly matches GeneLab RCP, so OSDR-provided STAR BAMs are reusable without re-alignment. |
| Ontology layer | GO (MF/BP/CC) + Plant Ontology (PO: anatomy/developmental) + Trait Ontology (TO) + MapMan bins. |
| FAIR deposit | Full package: GitHub repo with Zenodo DOI badge + Zenodo deposit with derived data tables, ML artifacts, ontology results, figures, RO-Crate metadata manifest. User uploads with their credentials. |
| Manuscript target | npj Microgravity (Nature portfolio). Full draft in Markdown + LaTeX source. |
| AS thresholds | FDR < 0.05 and |ΔPSI| > 0.1 for high-confidence events. |
| ML validation | Leave-one-study-out cross-validation (LOSOCV) as primary; nested k-fold as secondary comparison. |
| Visualization | Full suite: 12 figure types (splice graphs, PSI heatmaps, volcano plots, UMAP, event-type bars, ontology dot plots, MapMan maps, network graphs, circos plots, feature importance, sashimi plots, upset plots). |
| Execution | Full system built + PoC run on 6-8 studies. |
| Deposit prep | Complete package prepared locally; user pushes to GitHub + uploads to Zenodo. |
| Manuscript | Full prose draft in Markdown + LaTeX (Nature template) with all figures/tables. |

---

## 3. System Architecture (9 Modules)

### Module 1: OSDR Data Retrieval (`osdr_retrieval/`)
- **`osdr_client.py`**: Python client for the OSDR Biological Data API (BDAPI v2)
  - `enumerate_plant_studies()`: Query all 633 OSD datasets, filter by `organism` field + `study assay technology type` containing "RNA-Seq". Returns 54 Arabidopsis + 3 crop studies with metadata.
  - `get_study_metadata(osd_id)`: Fetch full ISA-Tab metadata (factors, tissue, genotype, condition, platform, mission).
  - `get_file_inventory(osd_id)`: List all files via BDAPI `/dataset/{id}/files/` endpoint, extract download URLs for `Aligned.sortedByCoord.out.bam` + `.bai` files.
  - `download_bams(osd_id, dest_dir)`: Download only the coordinate-sorted BAMs (not toTranscriptome BAMs) + index files. Stream to disk, verify MD5.
  - `parse_sample_metadata(osd_id)`: Extract per-sample factor values (Spaceflight/Ground Control, Light/Dark, tissue, genotype, replicate) from ISA-Tab metadata.
  - `build_contrast_matrix(osd_id)`: Auto-detect primary contrast (flight vs ground) and secondary contrasts (tissue, genotype, light) from factor types. Output rMATS-compatible b1.txt sample group definitions.
- **`study_catalog.json`**: Cached enumeration of all 57 plant RNA-seq studies with organism, title, factor types, sample counts, BAM availability, and download URLs.
- **`reference_files.py`**: Download and prepare Ensembl Plants release 48 references:
  - Arabidopsis: `Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz` + `Arabidopsis_thaliana.TAIR10.48.gtf.gz`
  - Brassica: `Brassica_rapa.Brapa_1.0.dna.toplevel.fa.gz` + `.Brapa_1.0.48.gtf.gz`
  - Brachypodium, tomato: corresponding Ensembl Plants release 48 files
  - Index GTF for rMATS; build SpliceGrapher gene model file.

### Module 2: Preprocessing & QC (`preprocessing/`)
- **`bam_qc.py`**: Parse STAR `Log.final.out` files (already in OSDR) for alignment rates, splice junction counts, mapping quality. Flag samples with <60% unique mapping.
- **`strandedness_check.py`**: Parse RSeQC `infer_experiment.py` results (already in OSDR MultiQC data) to confirm strandedness. Set rMATS `--libType` accordingly (unstranded/reversely/forwardly stranded).
- **`sample_sheet_builder.py`**: For each study, build rMATS b1.txt files defining replicate groups (b1 = treatment/flight, b2 = control/ground) from parsed OSDR metadata. Handle multi-factor studies by defining the primary contrast and noting secondary factors as covariates.
- **`qc_report.py`**: Compile per-study QC summary table (sample count, alignment rate, strandedness, read depth, replicate balance).

### Module 3: Differential Splicing (`diff_splicing/`)
- **`run_rmats.py`**: For each study + primary contrast:
  - Run `rmats.py --b1 flight_b1.txt --b2 ground_b1.txt --gtf Araport11.gtf --bam` (paired-end, strandedness from QC, `--statMethod` O, `--readLength` from STAR logs).
  - Produces 5 output files: `SE.MATS.JC.txt`, `A5SS.MATS.JC.txt`, `A3SS.MATS.JC.txt`, `MXE.MATS.JC.txt`, `RI.MATS.JC.txt`.
  - Each event: GeneID, geneSymbol, chr, strand, exon coordinates, IncLevelDifference (ΔPSI), FDR, inclusion/skipping counts per sample.
- **`rmats_postprocess.py`**:
  - Merge events across all studies into a unified event catalog (coordinate-based deduplication).
  - Apply thresholds: FDR < 0.05 AND |ΔPSI| > 0.1 → high-confidence differentially spliced events.
  - Build PSI matrix (events × samples) with NaN for missing events.
  - Annotate each event: gene ID, AGI code, event type, genomic coordinates, splice junction positions.
- **`nmd_prediction.py`**: For each AS event, predict whether the splicing change introduces a premature termination codon (PTC) → NMD susceptibility. Uses exon coordinates + ORF prediction. Classifies events as NMD-sensitive vs NMD-insensitive (functional consequence layer).
- **`domain_overlap.py`**: Map AS events to protein domain boundaries (InterPro/Pfam from Araport11 annotation). Classify: in-domain, domain-disrupting, domain-flanking, inter-domain. Feeds the functional impact ML model.

### Module 4: SpliceGrapher Visualization (`splice_grapher/`)
- **`setup_py2_env.sh`**: Create isolated Python 2.7 conda environment:
  - `conda create -n splicegrapher python=2.7`
  - Install SpliceGrapher from SourceForge, PyML 0.7.9, pysam 0.9, matplotlib 1.x
  - Download pre-built Arabidopsis splice site classifiers (SpliceGrapher ships these for >100 species including Arabidopsis).
- **`run_splicegrapher.py`** (runs in Py2 env):
  - `predict_graphs.py`: Generate splice graph predictions from BAM + Araport11 GTF for all genes with read coverage. Output `.gff` splice graph files.
  - `splicegrapher.py`: Render splice graph diagrams (exons as nodes, introns as edges, read depth overlay) for top differentially spliced loci.
  - Generate graphs for: (a) top 50 loci by |ΔPSI|, (b) top 20 loci by cross-study recurrence, (c) loci with known functional relevance from ontology enrichment.
- **`graph_selection.py`**: Select loci for visualization based on: statistical significance, effect size, cross-study recurrence, functional annotation (from ontology layer), and ML feature importance.

### Module 5: Machine Learning (`ml/`)
- **`feature_matrix.py`**: Build ML-ready feature matrix from PSI values:
  - Primary features: PSI values (events × samples), with k-NN imputation for missing values (k=5, distance = 1-|PSI correlation|).
  - Secondary features: event type (one-hot encoded), gene-level summaries (event count per gene, mean |ΔPSI| per gene).
  - Metadata labels: condition (flight/ground), study ID, tissue, genotype, environment type.
- **`classifier_flight_vs_ground.py`**:
  - Models: Random Forest, Gradient Boosting (XGBoost), linear SVM — compared.
  - **LOSOCV (primary)**: For each of N studies, train on N-1 studies, test on held-out study. Report per-study AUC, F1, accuracy. Aggregate: mean ± std across held-out studies.
  - **Nested k-fold (secondary)**: 5-fold outer, 3-fold inner, pooled across all samples. Report for comparison with LOSOCV.
  - Feature selection inside each CV fold (no leakage): variance filter + mutual information pre-filter to top 500 features, then model-specific importance.
  - Output: per-study held-out performance table, global feature importance ranking (SHAP values for tree models), top discriminating splicing events.
- **`clustering.py`**:
  - UMAP (n_neighbors=15, min_dist=0.1) on sample PSI profiles → 2D embedding, colored by condition/study/tissue/genotype.
  - Hierarchical clustering (Ward linkage, correlation distance) of genes by cross-study ΔPSI patterns → gene splicing modules.
  - Metrics: silhouette score, adjusted Rand index vs known labels (condition, study, tissue).
  - Output: UMAP scatter plots, dendrograms, cluster assignment tables, per-cluster enriched event types.
- **`functional_impact_model.py`**:
  - Target: predict whether an AS event disrupts a protein domain (binary, from `domain_overlap.py`).
  - Features: event type, exon length, position relative to domain boundaries, NMD prediction, ΔPSI magnitude, cross-study recurrence.
  - Validation: 5-fold CV with held-out loci (group CV by gene to prevent leakage of paralogous events).
  - Models: RF, GBM compared.
  - Output: per-event functional impact probability, feature importance, calibration plot.

### Module 6: Ontology Enrichment (`ontology/`)
- **`go_enrichment.py`**: GO over-representation analysis (ORA) + GSEA on differentially spliced gene sets (hypergeometric test, BH-FDR < 0.05). Uses `goatools` (Python) with Arabidopsis GO annotations from TAIR/GAF file. Three domains: MF, BP, CC.
- **`po_enrichment.py`**: Plant Ontology enrichment (anatomy terms: root, shoot, leaf, meristem; developmental stage terms). Uses PO OBO file + Arabidopsis PO annotations. Fisher's exact test, BH-FDR.
- **`to_enrichment.py`**: Trait Ontology enrichment (plant traits: stress tolerance, growth, development). Uses TO OBO file + trait annotations.
- **`mapman_mapping.py`**: MapMan bin assignment for differentially spliced genes (using Arabidopsis mapping file from MapMan). Pathway-level enrichment: which metabolic/regulatory pathways are enriched for splicing changes. Generates MapMan pathway maps with splicing-change overlay.
- **`cross_ontology.py`**: Integrate across GO + PO + TO + MapMan to build the adaptation narrative:
  - Identify genes/pathways consistently enriched across multiple ontologies.
  - Build a gene-to-function-to-adaptation mapping table.
  - Categorize splicing changes by biological process: stress response, cell wall remodeling, photosynthesis, hormone signaling, DNA repair, calcium signaling, oxidative stress.

### Module 7: Visualization Suite (`visualization/`)
All figures saved as SVG (primary) + PNG (secondary). Colorblind-friendly palettes. Liberation Sans font.

1. **Splice graph diagrams** (SpliceGrapher): Per-locus splice models for top differentially spliced genes, showing exon nodes, intron edges, read depth, and condition-specific splice junction usage.
2. **PSI heatmaps**: Clustered events × samples heatmap (rows = top differentially spliced events, columns = samples grouped by condition/study), color = PSI value.
3. **Volcano plots**: ΔPSI vs -log10(FDR) per study, with event-type coloring and labeled top loci.
4. **UMAP/t-SNE scatter**: Sample splicing profiles reduced to 2D, colored by condition (flight/ground), study, tissue, genotype. Shows whether splicing profiles cluster by environment.
5. **Event-type bar charts**: Stacked bar of SE/A5SS/A3SS/MXE/RI counts per study, showing which event types dominate in each environment.
6. **Ontology dot plots**: GO/PO/TO enrichment as dot plots (x = gene ratio, y = term, size = count, color = FDR).
7. **MapMan pathway maps**: Splicing changes overlaid on MapMan metabolic pathway diagrams, highlighting affected pathways.
8. **Network graphs**: Splicing-event → gene → function → adaptation network (Cytoscape-style layout via networkx), showing how splicing changes propagate to physiological processes.
9. **Circos plots**: Genomic distribution of differentially spliced events across Arabidopsis chromosomes, with tracks for event type, ΔPSI, and gene density.
10. **Feature importance plots**: Top discriminating splicing events from ML classifier (SHAP summary plot + bar chart of mean |SHAP|).
11. **Sashimi plots**: Read coverage + splice junction counts for top loci, flight vs ground overlay (via ggsashimi or matplotlib).
12. **Upset plots**: Overlap of differentially spliced genes across studies, showing core conserved vs environment-specific splicing changes.

### Module 8: FAIR Deposit (`fair_deposit/`)
- **`github_repo/`**: Complete code repository
  - All 7 analysis modules with docstrings and README per module
  - Top-level `README.md` with installation, usage, and pipeline diagram
  - `LICENSE` (MIT for code)
  - `environment.yml` (conda) + `requirements.txt` (pip) for Python 3 deps
  - `requirements_py2.txt` for SpliceGrapher Python 2.7 env
  - `CITATION.cff` (citation metadata)
  - `.github/workflows/ci.yml` (basic lint + import test)
  - Zenodo DOI badge in README (activates on first GitHub release)
  - `config/study_catalog.json` (all 57 plant RNA-seq studies)
  - `config/reference_config.yaml` (Ensembl Plants release 48 paths)
  - `snakemake/Snakefile` or `nextflow/main.nf` for full-pipeline orchestration (designed for all 54 studies, executable on cluster)
- **`zenodo_bundle/`**:
  - `derived_data/psi_matrix.csv` (events × samples PSI matrix)
  - `derived_data/diff_splicing_events.csv` (all significant events with annotations)
  - `derived_data/rmats_results/` (per-study rMATS output)
  - `derived_data/ml_model_artifacts/` (trained models, SHAP values, CV results)
  - `derived_data/ontology_results/` (GO/PO/TO/MapMan enrichment tables)
  - `derived_data/nmd_predictions.csv`
  - `derived_data/domain_overlap.csv`
  - `figures/` (all 12 figure types as SVG + PNG)
  - `tables/` (all manuscript tables as CSV)
  - `metadata/ro-crate-metadata.json` (RO-Crate 1.1 manifest linking all artifacts with DataCite metadata)
  - `metadata/datacite.json` (Zenodo deposit metadata: title, creators, description, keywords, license)
  - `README_deposit.md` (deposit description, file inventory, reuse instructions)
- **RO-Crate manifest**: JSON-LD following RO-Crate 1.1 specification, linking all digital objects with provenance, describing the workflow, inputs (OSDR study IDs), outputs, and software versions.

### Module 9: Manuscript (`manuscript/`)
- **`npj_microgravity_manuscript.md`**: Full prose draft in Markdown
  - Abstract (~250 words)
  - Introduction (AS in plant stress adaptation, spaceflight as complex stressor, OSDR resource, knowledge gap)
  - Results (system overview, AS landscape, splice graph models, ML classification, ML clustering, functional impact prediction, ontology enrichment, cross-species comparison)
  - Discussion (splicing as adaptation layer, key loci/pathways, ML insights, limitations, future directions)
  - Methods (detailed: data retrieval, rMATS, SpliceGrapher, ML, ontology, visualization, FAIR deposit)
  - References (inline citations, .bib file)
  - Figure captions (all 12 figures)
  - Table content (all tables)
- **`npj_microgravity.tex`**: LaTeX source using Nature portfolio template
  - `nature.cls` or `sn-jnl.cls` class file
  - Figures integrated with `\includegraphics`
  - Tables formatted with `booktabs`
  - Bibliography via `natbib` + `references.bib`
- **`references.bib`**: BibTeX bibliography
- **`figures/`**: All figures at 300 DPI (PNG) + vector (SVG/PDF) with figure numbers
- **`tables/`**: Supplementary tables as separate CSV files

---

## 4. PoC Study Selection (6-8 studies)

Selected to span the major environmental signal types in OSDR plant data:

| Study | Environment | Contrast | Samples | Rationale |
|---|---|---|---|---|
| OSD-120 | Spaceflight (ISS) | Flight vs Ground, root/shoot, Col-0/WS, light/dark | 36 | Rich multi-factor design, well-annotated |
| OSD-37 | Spaceflight (ISS) | Flight vs Ground, 4 ecotypes (Col-0, Ws-2, Ler-0, C24) | 56 | Cross-genotype splicing comparison |
| OSD-314 | Microgravity + Mars gravity | Microgravity vs Mars-g vs 1g, red light stimulation | 17 | Partial gravity splicing response |
| OSD-678 | Spaceflight (ISS) | Flight vs Ground, Col-0, light/dark | 36 | Clean flight-vs-ground with light factor |
| OSD-46 | Gamma radiation + HZE | Irradiated vs control | ~12 | Radiation stress splicing |
| OSD-134 | Low atmospheric pressure | Low pressure vs ambient | ~12 | Low-pressure stress splicing |
| OSD-251 | Fractional gravity | Partial-g vs 1g, blue light | ~12 | Partial gravity + light quality |
| OSD-59 | Spaceflight (ISS) | Flight vs Ground, Brassica/Mizuna | 4 | Cross-species case study (comparative) |

Total PoC: ~185 samples across 8 studies, spanning spaceflight, microgravity, partial gravity, radiation, low pressure, light quality, and cross-species.

---

## 5. Data Flow

```
OSDR BDAPI API
    │
    ▼
[Module 1: Retrieval] ──► study_catalog.json + BAM downloads + sample metadata
    │
    ▼
[Module 2: QC] ──► QC report + strandedness + sample sheets (b1.txt)
    │
    ▼
[Module 3: rMATS] ──► per-study AS events (SE/A5SS/A3SS/MXE/RI) + PSI matrix
    │                    + NMD predictions + domain overlap
    │
    ├──► [Module 4: SpliceGrapher] ──► splice graph diagrams (Py2 env)
    │
    ├──► [Module 5: ML] ──► classifier (LOSOCV) + clustering + functional impact
    │
    └──► [Module 6: Ontology] ──► GO + PO + TO + MapMan enrichment
                │
                ▼
         [Module 7: Visualization] ──► 12 figure types (SVG + PNG)
                │
                ▼
         [Module 8: FAIR Deposit] ──► GitHub repo + Zenodo bundle + RO-Crate
                │
                ▼
         [Module 9: Manuscript] ──► Markdown + LaTeX (npj Microgravity)
```

---

## 6. Compute/Resource Estimate

### PoC execution (6-8 studies, ~185 samples)

| Step | Input size | RAM | Disk | CPU | Est. runtime | Basis |
|---|---|---|---|---|---|---|
| BAM download (per study) | 1-2 GB/BAM × 12-56 BAMs | 1 GB | 20-80 GB peak | 1 | 20-60 min/study | OSDR download speed + Arabidopsis BAM size (small genome) |
| rMATS (per study) | 12-56 BAMs | 4-8 GB | 10-20 GB | 8 | 10-30 min/study | rMATS on Arabidopsis (~27k genes), junction-count mode |
| SpliceGrapher (50 loci) | 1 BAM + GTF | 2 GB | 5 GB | 2 | 60-120 min | Py2 predict_graphs + render, per-locus |
| ML (LOSOCV, 3 models) | ~10k features × 185 samples | 8 GB | 2 GB | 8 | 30-60 min | RF/GBM/SVM, 8 folds, feature selection per fold |
| Ontology enrichment | ~2k genes | 4 GB | 1 GB | 4 | 15-30 min | goatools + PO/TO Fisher tests + MapMan mapping |
| Visualization (12 types) | processed data | 4 GB | 5 GB | 4 | 30-60 min | matplotlib/seaborn/networkx/circos |
| **Total PoC** | — | **8-16 GB peak** | **~100 GB peak** | **8** | **~4-8 hours** | Sequential study processing |

### Execution strategy
- **worker-0 (default, 8 CPU / 32 GB)**: Retrieval, QC, ML, ontology, visualization, manuscript generation.
- **ManageMachine (right-sized, 8 CPU / 32 GB)**: rMATS processing — download one study's BAMs to /workspace, run rMATS, save results to /mnt/shared-workspace, delete BAMs, repeat. Keeps peak disk bounded.
- **SpliceGrapher**: Python 2.7 conda env on worker-0 (CPU-light, only visualization).
- **Checkpointing**: After each study's rMATS completes, save results to `/mnt/shared-workspace/checkpoints/`. If sandbox dies, resume from last completed study.
- **Full 54-study run**: Same code, orchestrated via Snakemake/Nextflow, designed for cluster execution. Not run in sandbox (would exceed 24h cap). Documented in manuscript Methods as "designed for and executable at full scale."

### Key risk mitigations
- **BAM download size**: Process one study at a time, delete BAMs after rMATS. Peak disk ~80 GB (one large study).
- **SpliceGrapher Py2 compatibility**: Test installation first in conda env. If PyML fails to install, fall back to `neograph` (Py2.7, simpler) or render splice graphs manually with matplotlib/networkx from rMATS junction coordinates.
- **rMATS with low replicate counts**: Some studies have 2 replicates/group. rMATS handles this but with reduced power. Flag in QC report; exclude from ML if insufficient.
- **Cross-study event merging**: Different studies may use slightly different STAR versions/parameters. Use coordinate-based event deduplication with ±2bp tolerance for junction positions.

---

## 7. Testing & Acceptance Criteria

- **OSDR retrieval**: Successfully enumerate all 57 plant RNA-seq studies; download + verify ≥1 BAM per study for PoC subset.
- **rMATS**: Produces non-empty output for all 5 event types for each PoC study; PSI matrix has <50% missing values.
- **SpliceGrapher**: Generates ≥20 readable splice graph diagrams for top loci. If Py2 env fails, fallback visualization produces equivalent diagrams.
- **ML classifier**: LOSOCV completes for all PoC studies; reports per-study AUC; feature importance identifies ≥10 named splicing events.
- **Clustering**: UMAP produces 2D embedding with visible condition separation (if signal exists); silhouette score reported.
- **Functional impact**: Model trains and produces per-event probabilities; calibration plot generated.
- **Ontology**: At least one GO/PO/TO term enriched at FDR < 0.05 across the PoC gene set; MapMan mapping covers >70% of differentially spliced genes.
- **Visualization**: All 12 figure types generated, non-blank, readable (verified via media output check).
- **FAIR deposit**: GitHub repo structure complete with all code + configs; Zenodo bundle contains all derived data + figures + RO-Crate manifest; metadata valid against DataCite schema.
- **Manuscript**: Full Markdown draft with all sections, ≥8 figures referenced, ≥4 tables; LaTeX compiles (or is structurally valid for manual compilation).

---

## 8. Assumptions

1. OSDR BDAPI remains accessible and returns the same study/file structure observed during planning.
2. Pre-aligned STAR BAMs in OSDR are aligned to Ensembl Plants release 48 TAIR10 + Araport11 (confirmed for GL-DPPD-7101-E pipeline; older studies may use earlier RCP versions with different references — will verify per-study and re-align from FASTQ if reference mismatch detected).
3. rMATS installs cleanly via conda/pip in the sandbox Python 3 environment.
4. SpliceGrapher installs in a Python 2.7 conda environment; if not, fallback to manual splice graph rendering from rMATS junction data.
5. PoC studies have ≥2 biological replicates per condition (will verify from metadata; exclude if not).
6. Arabidopsis GO/PO/TO/MapMan annotation files are freely available from TAIR/PlantOntology/MapMan.
7. User has GitHub and Zenodo accounts and will upload the prepared deposit package.
8. npj Microgravity accepts LaTeX submissions via the Nature portfolio submission system.
