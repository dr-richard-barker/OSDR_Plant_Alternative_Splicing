# Plant Spaceflight Alternative Splicing Analysis — FAIR Zenodo Deposit Guide

This directory contains the complete metadata, data dictionaries, and FAIR assets required for depositing our alternative splicing pipeline results and models to **Zenodo** and creating a peer-reviewed **RO-Crate** manifest.

## Deposit Packaged Contents

- `CITATION.cff` — Citation File Format metadata, detailing preferred citation formats, contributors, and repository coordinates.
- `LICENSE` — Software code and documentation are released under the open-source **MIT License**.
- `zenodo.json` — Pre-reserved publication schema for the Zenodo REST API, including metadata keywords, author affiliations (NASA GeneLab / OSDR), description, and target community (`nasa-genelab`).
- `ro-crate-metadata.json` — A lightweight, compliant **RO-Crate (v1.1)** metadata manifest linking all physical files in the repository (such as the pooled feature matrices, individual rMATS outputs, UMAP projections, GO over-representation tables, and splice graphs) with formal semantic schemas (W3C/DataCite/Schema.org).
- `data_dictionary.json` — Fully annotated schemas mapping every metric, statistics column, and PSI matrix coordinate to precise mathematical and biological definitions.

---

## Detailed Data Dictionary & Schema Mapping

The following specifications detail the core tables included in the Zenodo deposit and available in our repository:

### 1. Splicing Feature Matrix (`multi_study_analysis_v4/feature_matrix.npy` & `feature_genes.tsv`)
- **feature_matrix.npy:** A binary float32 numpy array of dimension `199 × 9,875` representing the imputed **Percent Spliced In (PSI)** values.
- **feature_genes.tsv:** A single-column text table mapping the 9,875 columns of the feature matrix to unique Arabidopsis locus identifiers (AGI codes).

### 2. Sample Metadata (`multi_study_analysis_v4/sample_metadata.tsv`)
- **Sample_ID:** Unique coordinate name of the sequencing library.
- **OSD_ID:** Open Science Data Repository accession (e.g., `OSD-314`, `OSD-120`).
- **Condition:** Categorical environment classification (`Treatment` for flight/altered gravity/irradiated/regolith vs `Control` for ground/1g/unirradiated/JSC-1A simulant).
- **Tissue:** Sample tissue origin (e.g., `Root`, `Shoot`).
- **Genotype:** Genetic background (e.g., `Col-0`, `Ws-2`, `phyD`).

### 3. Significant Splicing Event Catalogs (`osd*_analysis/OSD-*_significant_events.tsv`)
- **GeneID:** Araport11 gene identifier.
- **geneSymbol:** Standard gene nomenclature.
- **chr / strand / start / end:** Genomic coordinates of the alternative splicing event.
- **FDR:** Benjamini-Hochberg false discovery rate from the rMATS-turbo likelihood ratio test.
- **IncLevelDifference (ΔPSI):** Splicing change between Treatment and Control. Positive values denote increased splice inclusion in space stress, negative denotes decreased inclusion.

### 4. Nonsense-Mediated Decay Heuristics (`osd*_analysis/OSD-*_nmd_predictions.tsv`)
- **Event_ID:** Coordinate-based unique event key.
- **Frameshift_Heuristic:** Divisibility check of the alternative sequence length by 3.
- **NMD_Sensitive:** Boolean field representing susceptibility to premature termination codons (PTC) leading to transcript decay.

---

## Instructions for Uploading to Zenodo

To complete the peer-reviewed deposition of these FAIR assets:

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
