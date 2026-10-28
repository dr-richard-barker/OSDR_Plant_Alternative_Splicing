# Plant Spaceflight Alternative Splicing Analysis — FAIR Zenodo Deposit Guide

This directory contains the metadata, data dictionary, and FAIR assets for depositing the alternative splicing results to **Zenodo**, with an **RO-Crate** manifest. TODO: the analysis code is not yet in the repository, and nothing has been deposited yet.

## Deposit Packaged Contents

- `CITATION.cff` — Citation File Format metadata, detailing preferred citation formats, contributors, and repository coordinates.
- `LICENSE` — Software code and documentation are released under the open-source **MIT License**.
- `zenodo.json` — Pre-reserved publication schema for the Zenodo REST API, including metadata keywords, author affiliations (NASA GeneLab / OSDR), description, and target community (`nasa-genelab`).
- `ro-crate-metadata.json` — **RO-Crate (v1.1)** metadata manifest listing the result files in the repository (pooled feature matrices, individual rMATS outputs, UMAP projections, GO over-representation tables, splice graphs and manuscript drafts).
- `data_dictionary.json` — Per-file columns, row counts, shapes and dtypes, generated from the files themselves by `tools/build_data_dictionary.py`, plus definitions of the main columns.

---

## Detailed Data Dictionary & Schema Mapping

The following specifications detail the core tables included in the Zenodo deposit and available in our repository:

### 1. Splicing Feature Matrix (`multi_study_analysis_v4/feature_matrix.npy` & `feature_genes.tsv`)
- **feature_matrix.npy:** A float64 NumPy array of shape `199 × 9,875` holding gene-level mean **Percent Spliced In (PSI)** values (mean-imputed).
- **feature_genes.tsv:** A single-column text table mapping the 9,875 columns of the feature matrix to unique Arabidopsis locus identifiers (AGI codes).

### 2. Sample Metadata (`multi_study_analysis_v4/sample_metadata.tsv`)
Columns: `sample_name`, `study_id`, `condition`, `is_flight` (199 rows).
- **sample_name:** Sample name as in the GeneLab file names.
- **study_id:** OSDR accession (e.g., `OSD-314`, `OSD-120`).
- **condition:** Condition label used for the contrast.
- **is_flight:** 1 = treatment (flight, altered gravity, irradiated or regolith), 0 = control (ground, 1 g, non-irradiated or JSC-1A simulant).
- Tissue and genotype are not in this table; per-study tissue (from OSDR) is listed in the top-level README.

### 3. Significant Splicing Event Catalogs (`osd*_analysis/OSD-*_significant_events.tsv`)
- **GeneID:** Araport11/TAIR10 gene identifier (AGI code).
- **geneSymbol:** Standard gene nomenclature.
- **chr / strand / start / end:** Genomic coordinates of the alternative splicing event.
- **FDR:** Benjamini-Hochberg false discovery rate from the rMATS-turbo likelihood ratio test.
- **IncLevelDifference (ΔPSI):** rMATS mean(IncLevel1) − mean(IncLevel2), i.e. sample group 1 minus group 2. Where group sizes differ (OSD-314, OSD-658, OSD-476, OSD-251), group 1 is the treatment group.

### 4. Nonsense-Mediated Decay Heuristics (`osd*_analysis/OSD-*_nmd_predictions.tsv`)
Two layouts exist. OSD-314 is tab-separated and OSD-120, OSD-37 and OSD-678 are comma-separated, all with `event_id, event_type, nmd_prediction, nmd_reason, frameshift`. OSD-658, OSD-476 and OSD-251 are tab-separated with `event_id, GeneID, event_type, NMD_prediction`.
- **nmd_prediction / NMD_prediction:** Heuristic class, not a measurement. SE/A5SS/A3SS/MXE events whose length is not divisible by 3 are classed NMD-sensitive, and every RI event is classed NMD-sensitive.
- **nmd_reason / frameshift:** The rule that produced the class, and whether a frameshift was inferred.

---

## Instructions for Uploading to Zenodo

To deposit these assets:

1. **Verify metadata formats:** Ensure `zenodo.json` has correct author affiliations and keyword schemas.
2. **Compress analysis assets:** Create a unified archive containing `docs/`, `multi_study_analysis_v4/`, the `osd*_analysis/` directories, and the contents of this `fair_deposit/` folder.
3. **Execute REST API Upload:**
   ```bash
   # Example shell command to upload through the Zenodo API using pre-reserved credentials
   curl -H "Content-Type: application/json" \
        -H "Authorization: Bearer $ZENODO_TOKEN" \
        -X POST \
        -d @fair_deposit/zenodo.json \
        https://zenodo.org/api/deposit/depositions
   ```
4. **Publish & Link DOI:** Associate the newly generated Zenodo DOI with the GitHub Pages repository README badge.
