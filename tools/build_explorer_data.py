#!/usr/bin/env python3
"""Regenerate the splice-explorer gene table (the STUDIES block) in docs/index.html.

Every per-gene field shown in the explorer is derived from a file, not typed by hand:

* symbol       -- Ensembl Plants display_name (tools/ensembl_gene_lookup.json);
                  "—" when Ensembl has no symbol beyond the AGI code.
* description  -- the ARAPORT_DETAILS block already in docs/index.html (checked here
                  against the Ensembl lookup: description, isoform count, canonical
                  isoform must all agree), followed by the gene's significant events
                  in that study (osd*_analysis/OSD-*_significant_events.tsv) and its
                  recurrence (multi_study_analysis_v4/recurrent_genes.tsv).
* recurrent    -- True iff the gene is listed in recurrent_genes.tsv.

The gene list and figure filenames per study are taken from the existing STUDIES
block, because they correspond to the splice-graph figures in docs/*/figures/.
Study labels (stressor, tissue) and sample counts are from the OSDR study records
(https://visualization.osdr.nasa.gov/biodata/api/v2/), read 2026-10-28.

Usage:  python3 tools/build_explorer_data.py          (rewrites docs/index.html)
"""
import csv
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
HTML = ROOT / "docs" / "index.html"
BEGIN = "// BEGIN GENERATED STUDIES (tools/build_explorer_data.py)"
END = "// END GENERATED STUDIES"

# From OSDR sample records (factor values / material type / organism part).
STUDY_LABELS = {
    "osd314": "Altered gravity on the ISS (µg and 0.34 g) vs 1 g on Earth · whole seedlings",
    "osd120": "ISS spaceflight vs ground (CARA flight) · root tips",
    "osd37": "ISS spaceflight vs ground, 4 ecotypes · whole seedlings",
    "osd251": "Fractional gravity on the ISS (µg–0.57 g vs in-flight 1 g), blue light · seedlings",
    "osd476": "Apollo 11/12/17 regolith vs JSC-1A simulant (ground) · aboveground tissue",
    "osd658": "Simulated GCR, 40 and 80 cGy vs non-irradiated (ground) · whole seedlings",
    "osd678": "ISS spaceflight vs ground (same CARA flight as OSD-120) · leaves",
}


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fmt_dpsi(x):
    return f"{float(x):+.3f}"


def fmt_fdr(x):
    x = float(x)
    return f"{x:.1e}" if x < 0.001 else f"{x:.3f}"


def main():
    html = HTML.read_text()
    ensembl = json.loads((ROOT / "tools" / "ensembl_gene_lookup.json").read_text())["genes"]

    m = re.search(r"const ARAPORT_DETAILS = (\{.*?\n\});", html, re.S)
    araport = json.loads(m.group(1))
    for gid, d in araport.items():
        e = ensembl[gid]
        assert d["description"] == e["description"], (gid, d["description"], e["description"])
        assert d["isoforms"] == e["n_transcripts"], (gid, d["isoforms"], e["n_transcripts"])
        assert d["canonical"] == e["canonical_transcript"], (gid, d["canonical"], e["canonical_transcript"])

    if BEGIN in html:
        start = html.index(BEGIN)
        stop = html.index(END) + len(END)
        block = html[start:stop]
    else:
        start = html.index("const STUDIES = [")
        stop = html.index("\n];", start) + len("\n];")
        block = html[start:stop]

    studies = []
    if BEGIN in html:  # previously generated block: plain JSON
        data = json.loads(block[block.index("const STUDIES = ") + len("const STUDIES = "):block.rindex(";")])
        for s in data:
            studies.append((s["id"], s["name"], s["samples"], [(g["id"], g["filename"]) for g in s["genes"]]))
    else:  # original hand-written block
        for sm in re.finditer(r'id: "(osd\d+)",\s*name: "([^"]+)",.*?samples: (\d+),\s*genes: \[(.*?)\n    \]', block, re.S):
            sid, name, samples, genes_block = sm.groups()
            genes = re.findall(r'\{ id: "([^"]+)",.*?filename: "([^"]+)"', genes_block)
            studies.append((sid, name, int(samples), genes))
    assert len(studies) == 7, len(studies)

    recurrent = {r["gene"]: r for r in read_tsv(ROOT / "multi_study_analysis_v4" / "recurrent_genes.tsv")}

    out_studies = []
    for sid, name, samples, genes in studies:
        num = sid[3:]
        sig = read_tsv(ROOT / f"{sid}_analysis" / f"OSD-{num}_significant_events.tsv")
        n1 = len(sig[0]["IJC_SAMPLE_1"].split(","))
        n2 = len(sig[0]["IJC_SAMPLE_2"].split(","))
        assert n1 + n2 == samples, (sid, n1, n2, samples)
        out_genes = []
        for gid, filename in genes:
            if gid.startswith("chr"):
                ev = [r for r in sig if gid.split(":")[1].split("-")[0] in r["event_id"]]
                assert len(ev) == 1, (gid, len(ev))
                agi = ev[0]["GeneID"]
                symbol = ensembl[agi]["display_name"] if agi in ensembl else agi
                desc = (f"Splice graph for the {ev[0]['event_type']} event {gid} in {agi}"
                        f" ({symbol}): ΔPSI {fmt_dpsi(ev[0]['IncLevelDifference'])}, FDR {fmt_fdr(ev[0]['FDR'])}.")
                out_genes.append({"id": gid, "symbol": f"{symbol} RI locus", "filename": filename,
                                  "description": desc, "recurrent": agi in recurrent})
                continue
            e = ensembl[gid]
            symbol = e["display_name"] if e["display_name"] and e["display_name"] != gid else "—"
            ev = [r for r in sig if r["GeneID"] == gid]
            assert ev, (sid, gid, "not significant in this study")
            evtxt = "; ".join(f"{r['event_type']} ΔPSI {fmt_dpsi(r['IncLevelDifference'])} (FDR {fmt_fdr(r['FDR'])})"
                              for r in sorted(ev, key=lambda r: float(r["FDR"])))
            desc = f"Araport11: {araport[gid]['description']}. Significant in {name}: {evtxt}."
            if gid in recurrent:
                rr = recurrent[gid]
                desc += f" Significant in {rr['n_studies']} of 7 studies ({rr['studies'].replace(',', ', ')})."
            out_genes.append({"id": gid, "symbol": symbol, "filename": filename,
                              "description": desc, "recurrent": gid in recurrent})
        out_studies.append({"id": sid, "name": name, "stressor": STUDY_LABELS[sid],
                            "samples": samples, "genes": out_genes})

    js = (BEGIN + "\n// Do not edit by hand; rerun the script.\nconst STUDIES = "
          + json.dumps(out_studies, indent=2, ensure_ascii=False) + ";\n" + END)
    HTML.write_text(html[:start] + js + html[stop:])
    n = sum(len(s["genes"]) for s in out_studies)
    print(f"wrote {len(out_studies)} studies, {n} gene entries, "
          f"{sum(g['recurrent'] for s in out_studies for g in s['genes'])} flagged recurrent")


if __name__ == "__main__":
    sys.exit(main())
