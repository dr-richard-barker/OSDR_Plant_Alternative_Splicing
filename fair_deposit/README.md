# Plant Spaceflight Alternative Splicing Analysis

## Overview

Differential alternative splicing (AS) analysis of plant RNA-seq studies from the
NASA Open Science Data Repository (OSDR), using rMATS-turbo for event quantification,
machine learning for pattern discovery, and multi-ontology functional annotation.

## Contents

- `code/` - Complete analysis pipeline (Python)
  - `osdr_retrieval/` - OSDR BDAPI client for study enumeration and BAM download
  - `preprocessing/` - QC, strandedness detection, rMATS sample sheet builder
  - `diff_splicing/` - rMATS runner, event merging, PSI matrix, NMD/domain prediction
  - `ml/` - LOSOCV classifier, UMAP clustering, functional impact prediction
  - `ontology/` - GO, Plant Ontology, Trait Ontology, MapMan enrichment
  - `visualization/` - 12-figure publication suite
  - `fair_deposit/` - Zenodo/RO-Crate package builder
- `results/` - Analysis outputs
  - `derived_data/` - Event catalogs, PSI matrices, ML results, enrichment tables
  - `figures/` - Publication figures (SVG + PNG)
- `config/` - Study catalog and PoC configuration

## Studies Analyzed (proof-of-concept)

### Completed rMATS + Full Analysis

| OSD ID | Organism | Environment | Contrast | Samples | Total Events | Sig Events |
|--------|----------|-------------|----------|---------|-------------|------------|
| OSD-314 | Arabidopsis | Microgravity/Mars-g | altered gravity vs 1g | 17 | 10,532 | 161 |
| OSD-120 | Arabidopsis | ISS spaceflight | flight vs ground | 36 | 13,570 | 32 |
| OSD-37 | Arabidopsis | ISS spaceflight (4 accessions) | flight vs ground | 56 | 20,590 | 10 |
| OSD-59 | Brassica rapa | ISS spaceflight | flight vs ground | 2 | 1,877 | 77 |
| OSD-678 | Arabidopsis | ISS spaceflight (3 genotypes × 2 light) | flight vs ground | 36 | 15,255 | 190 |
| OSD-658 | Arabidopsis | Simulated GCR radiation | irradiated vs control | 14 | 9,932 | 43 |
| OSD-476 | Arabidopsis | Lunar regolith | Apollo regolith vs JSC-1A | 20 | 14,766 | 306 |
| OSD-251 | Arabidopsis | Fractional gravity (blue light) | altered gravity vs 1g | 20 | 15,220 | 86 |

### Multi-Study Analysis (7 Arabidopsis studies)

| Metric | Value |
|--------|-------|
| Total samples | 199 (112 treatment, 87 control) |
| Feature matrix | 199 × 9,875 genes |
| LOSOCV (7-study) | RF AUC=0.412±0.112, GBM AUC=0.412±0.133, XGBoost AUC=0.412±0.112 |
| Nested 5-fold CV | RF AUC=0.598±0.084 |
| UMAP ARI (study) | 0.995 (batch effect dominates) |
| UMAP ARI (treatment) | -0.002 |
| Cross-study recurrence | 93 genes in ≥2 studies; AT4G17310 in 4 studies; 8 genes in 3 studies |
| SHAP top predictors | AT1G56220, AT5G01770, AT1G48360, AT3G15620, AT3G62190 |

## Key Results

- **RI dominance**: Retained intron events constitute 57-67% of significant events across all 7 Arabidopsis studies
- **NMD sensitivity**: 79-86% of significant events are NMD-sensitive across all 7 Arabidopsis studies (regulatory role)
- **Circadian enrichment**: Circadian rhythm regulation enriched in ALL 4 ISS flight studies (OSD-314 FDR=0.032, OSD-120 FDR=0.097, OSD-37 via RVE1, OSD-678 FDR=0.008 FE=18.6x) — most robust finding within flight-study subset
- **Stressor-specific GO profiles**: GCR radiation (OSD-658) → phospholipid metabolism; lunar regolith (OSD-476) → photosynthesis/cell cycle; fractional gravity (OSD-251) → superoxide response (FDR=0.001, FE=488x)
- **Cross-study recurrence**: 93 genes significant in ≥2 studies; AT4G17310 in 4 studies (OSD-251, OSD-314, OSD-476, OSD-658) — most robust splicing target across gravity, radiation, and regolith; 8 genes in 3 studies (SOS4, CRK28, QWRF3, etc.)
- **Largest pairwise overlaps**: OSD-678 ∩ OSD-476 = 26 genes; OSD-314 ∩ OSD-476 = 22 genes; OSD-314 ∩ OSD-678 = 18 genes
- **Effect-size threshold importance**: OSD-37 (56 samples) detected 175 FDR-significant events but only 10 with |ΔPSI|>0.1 (median |ΔPSI|=0.025)
- **ML classification**: Nested 5-fold CV AUC=0.598 (7-study pooled); 7-study LOSOCV near-random (RF/GBM/XGBoost=0.412)
- **Batch effect dominance**: UMAP shows study identity almost perfectly separates samples (ARI=0.995) while treatment/control does not (ARI=-0.002)
- **SHAP biomarkers**: AT1G56220, AT5G01770, AT1G48360, AT3G15620, AT3G62190 are top 7-study splicing-based predictors
- **OSD-476 most splicing-responsive**: 306 significant events in 257 genes (largest of any study)

## Methods

- **Differential splicing**: rMATS-turbo v4.3.0, FDR < 0.05, |dPSI| > 0.1
- **Event types**: SE, A5SS, A3SS, MXE, RI
- **Reference**: Ensembl Plants release 48 (TAIR10 for Arabidopsis, Brapa_1.0 for Brassica)
- **ML**: LOSOCV (RF/GBM/XGBoost), UMAP clustering, SHAP feature importance
- **Ontology**: GO (MF/BP/CC), Plant Ontology, Trait Ontology, MapMan bins
- **Functional impact**: NMD susceptibility (frameshift heuristic), protein domain overlap

## Requirements

- Python 3.10+
- rMATS-turbo (conda: `conda create -n rmats -c bioconda -c conda-forge rmats python=3.10`)
- Key packages: pandas, numpy, scipy, scikit-learn, xgboost, shap, umap-learn, pysam, seaborn, matplotlib

## Usage

```bash
# 1. Enumerate OSDR plant studies
python code/osdr_retrieval/osdr_client.py

# 2. Download BAMs for a study
python -c "from code.osdr_retrieval.osdr_client import download_bams; download_bams('OSD-314', 'data/OSD-314')"

# 3. Run rMATS
python code/diff_splicing/run_rmats.py --poc-config config/poc_studies.json

# 4. Post-process events
python code/diff_splicing/rmats_postprocess.py

# 5. Run ML analysis
python code/ml/analysis.py

# 6. Generate figures
python code/visualization/generate_figures.py
```

## License

MIT License - see LICENSE file

## Citation

See CITATION.cff

## Data Source

All RNA-seq data was obtained from the NASA Open Science Data Repository (OSDR):
https://osdr.nasa.gov/biodata/api/
