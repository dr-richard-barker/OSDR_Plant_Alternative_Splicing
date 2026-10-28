#!/usr/bin/env python3
"""Regenerate fair_deposit/data_dictionary.json from the data files actually in the repository.

Shapes, dtypes, column names and row counts are read from the files; per-file descriptions
are templates filled with counts computed here. Usage: python3 tools/build_data_dictionary.py
"""
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "fair_deposit" / "data_dictionary.json"

# Files the previous dictionary described but which are not in this repository.
NOT_IN_REPO = ["event_catalog.csv", "significant_events.csv", "psi_matrix.csv", "study_summary.csv",
               "nmd_predictions.csv", "domain_overlap.csv", "go_enrichment.csv", "po_enrichment.csv",
               "to_enrichment.csv", "mapman_enrichment.csv", "ml_results/losocv_rf_per_study.csv",
               "ml_results/losocv_rf_top_features.csv", "ml_results/umap_clustering.csv"]

COLUMN_DEFINITIONS = {
    "GeneID": "Araport11/TAIR10 gene identifier (AGI code) from the Ensembl Plants annotation",
    "geneSymbol": "Gene symbol from the annotation (empty when none)",
    "chr / strand": "Chromosome and strand of the event",
    "IJC_SAMPLE_1 / SJC_SAMPLE_1 / IJC_SAMPLE_2 / SJC_SAMPLE_2": "rMATS inclusion / skipping junction counts per replicate (comma-separated) for sample groups 1 and 2",
    "PValue / FDR": "rMATS-turbo likelihood-ratio p-value and Benjamini-Hochberg FDR",
    "IncLevel1 / IncLevel2": "Per-replicate PSI for sample groups 1 and 2 (comma-separated)",
    "IncLevelDifference": "rMATS mean(IncLevel1) - mean(IncLevel2); group 1 is the treatment group where group sizes differ (OSD-314, OSD-658, OSD-476, OSD-251)",
    "dPSI": "Copy of IncLevelDifference",
    "event_type": "SE, A5SS, A3SS, MXE or RI",
    "event_id": "Coordinate-based event key: <type>_<chr>_<strand>_<coordinates>",
    "nmd_prediction / NMD_prediction": "Heuristic class, not a measurement: SE/A5SS/A3SS/MXE length not divisible by 3 -> NMD-sensitive; every RI -> NMD-sensitive",
    "nmd_reason / frameshift": "Rule that produced the class / whether a frameshift was inferred",
    "go_id / GO_ID, go_name / term, namespace / ontology": "Gene Ontology term and namespace",
    "n_genes_test / gene_count, genes": "Number and list of significant genes annotated to the term",
    "fold_enrichment, pvalue / p_value, fdr_bh / FDR": "Fisher's exact test over-representation statistics; FDR is Benjamini-Hochberg over the terms listed in the file",
    "is_flight": "1 = treatment (flight, altered gravity, radiation or regolith), 0 = control",
}


def tracked():
    out = subprocess.check_output(["git", "ls-files"], cwd=ROOT).decode().split("\n")
    return [p for p in out if p.endswith((".tsv", ".npy")) and (p.startswith("osd") or p.startswith("multi_study"))]


def tsv_entry(path):
    df = pd.read_csv(ROOT / path, sep=None, engine="python") if path.endswith("nmd_predictions.tsv") else \
        pd.read_csv(ROOT / path, sep="\t", low_memory=False)
    e = {"format": "tab-separated text" if "\t" in (ROOT / path).read_text().split("\n", 1)[0] else "comma-separated text (despite .tsv extension)",
         "n_rows": int(len(df)), "columns": list(df.columns)}
    name = Path(path).name
    study = name.split("_")[0] if name.startswith("OSD-") else None
    if name.endswith("_all_events.tsv"):
        e["description"] = f"All {len(df):,} rMATS-turbo events for {study}; by type: " + \
            ", ".join(f"{k} {v:,}" for k, v in df.event_type.value_counts().items())
    elif name.endswith("_significant_events.tsv"):
        e["description"] = f"{len(df)} significant {study} events (FDR < 0.05 and |IncLevelDifference| > 0.1); by type: " + \
            ", ".join(f"{k} {v}" for k, v in df.event_type.value_counts().items())
    elif name.endswith("_nmd_predictions.tsv"):
        col = "nmd_prediction" if "nmd_prediction" in df else "NMD_prediction"
        k = int((df[col].str.lower().str.replace("-", "_") == "nmd_sensitive").sum())
        e["description"] = f"Heuristic NMD class for the {len(df)} significant {study} events ({k} classified NMD-sensitive); see column_definitions"
    elif name.endswith("_go_enrichment.tsv"):
        p = "pvalue" if "pvalue" in df else "p_value"
        f = "fdr_bh" if "fdr_bh" in df else "FDR"
        e["description"] = f"GO over-representation results for {study}: {len(df)} terms listed, {(df[p] < 0.05).sum()} with p < 0.05, {(df[f] < 0.25).sum()} with FDR < 0.25"
    elif name.endswith("_mapman_enrichment.tsv"):
        e["description"] = f"MapMan-bin over-representation results for {study} (GO-derived bin mapping): {len(df)} bins listed"
    elif name.endswith("_psi_matrix.tsv"):
        e["description"] = f"PSI matrix for {study}: {len(df):,} events x {len(df.columns) - 3} samples"
    elif name.endswith("_umap_clustering.tsv"):
        e["description"] = f"UMAP coordinates and cluster labels for {study} ({len(df)} samples)"
    elif name == "recurrent_genes.tsv":
        vc = df.n_studies.value_counts().sort_index()
        e["description"] = f"{len(df)} genes significant in >= 2 of 7 studies (" + ", ".join(f"{v} in {k}" for k, v in vc.items()) + ")"
    elif name == "cross_study_recurrence.tsv":
        e["description"] = "Pairwise overlap of significant genes between studies"
    elif name == "losocv_results.tsv":
        m = df.groupby("model").auc.agg(["mean", "std"])
        e["description"] = "Leave-one-study-out classification results; mean AUC +/- SD: " + \
            ", ".join(f"{k} {r['mean']:.3f} +/- {r['std']:.3f}" for k, r in m.iterrows())
    elif name == "nested_kfold_results.tsv":
        e["description"] = f"Nested 5-fold CV (Random Forest) AUC per fold; mean {df.auc.mean():.3f} +/- {df.auc.std():.3f}"
    elif name == "umap_clustering.tsv":
        e["description"] = f"UMAP coordinates and cluster labels ({len(df)} samples)"
    elif name == "shap_importance.tsv":
        e["description"] = f"Mean |SHAP| per gene ({len(df)} genes), sorted descending"
    elif name == "sample_metadata.tsv":
        e["description"] = f"Sample metadata ({len(df)} samples; " + ", ".join(f"{k} {v}" for k, v in df.study_id.value_counts().items()) + ")"
    elif name == "feature_genes.tsv":
        e["description"] = f"Gene identifiers for the {len(df):,} columns of feature_matrix.npy"
    else:
        e["description"] = "TODO: describe"
    return e


def npy_entry(path):
    a = np.load(ROOT / path, allow_pickle=True)
    e = {"format": "NumPy .npy" + (" (object array; load with allow_pickle=True)" if a.dtype == object else ""),
         "dtype": str(a.dtype), "shape": list(a.shape)}
    name = Path(path).name
    if name == "feature_matrix.npy":
        e["description"] = f"Gene-level PSI feature matrix, {a.shape[0]} samples x {a.shape[1]:,} genes"
    elif name == "labels_flight.npy":
        e["description"] = f"Binary treatment labels ({int((a == 1).sum())} treatment, {int((a == 0).sum())} control)"
    elif name == "study_ids.npy":
        e["description"] = "Study accession per sample: " + ", ".join(sorted(set(map(str, a))))
    return e


def main():
    files = {}
    for p in sorted(tracked()):
        files[p] = npy_entry(p) if p.endswith(".npy") else tsv_entry(p)
        if p.startswith("multi_study_analysis/"):
            files[p]["description"] = "Earlier 2-study run (OSD-314 + OSD-120), superseded by multi_study_analysis_v4. " + files[p]["description"]
    out = {
        "description": "Data dictionary for the alternative-splicing result files in this repository, generated by tools/build_data_dictionary.py from the files themselves.",
        "column_definitions": COLUMN_DEFINITIONS,
        "files": files,
        "not_in_repository": {"files": NOT_IN_REPO,
                              "note": "TODO: described in an earlier version of this dictionary but not present in the repository."},
    }
    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n")
    print(f"wrote {len(files)} file entries to {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
