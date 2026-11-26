#!/usr/bin/env python3
"""R1: run-level summary across all samples of a run -> R1_run_qc.summary.json + run_qc.pdf
(rendered with the same step_report layout: plain-language page, then technical pages)."""
import argparse, collections, datetime, glob, json, os
import yaml
from step_report import build

STEPS = ["00_input", "01_read_qc", "02_species", "03_alignment", "04_coverage", "05_variants", "06_deletions",
         "07_annotation", "08_lineage", "09_resistance", "10_mixed_infection", "11_classification", "12_report"]


def main():
    ap = argparse.ArgumentParser()
    for k in ("summaries", "text", "version", "out_pdf"):
        ap.add_argument("--" + k, required=True)
    a = ap.parse_args()
    by = collections.defaultdict(dict)
    run = None
    for p in glob.glob(os.path.join(a.summaries, "*.json")):
        d = json.load(open(p))
        by[d["sample_id"]][d["step"]] = d
        run = d["run_id"]
    rows, status_rows, n_rep = [], [], 0
    for sid in sorted(by):
        st = by[sid]
        rep = st.get("12_report")
        last = max((x for x in STEPS if x in st), key=STEPS.index)
        if rep:
            v = rep["values"]
            qc = v["qc"]
            n_rep += qc != "FAIL"
            res = ", ".join(v["detected"]) or "none"
            rows.append([sid, qc, v.get("category") or "-", res, ", ".join(v["indeterminate"]) or "-"])
        else:
            rows.append([sid, "FAIL", f"stopped at {last}", "-", "-"])
        status_rows.append([sid] + [{"PASS": "P", "WARN": "C", "FAIL": "F"}.get(st[x]["status"], "-")
                                    if x in st else "-" for x in STEPS[:12]])
    n = len(rows)
    summ = {"schema": "mtb-step-summary/1", "sample_id": f"run {run}", "run_id": run, "step": "R1_run_qc",
            "status": "PASS" if n_rep == n else "WARN", "pipeline_version": a.version,
            "generated": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
            "values": {"n": n, "n_rep": n_rep, "n_fail": n - n_rep},
            "metrics": [{"name": "Samples with a report", "plain": "Samples with a usable result",
                         "value": f"{n_rep} of {n}", "threshold": "all", "status": "PASS" if n_rep == n else "WARN"}],
            "messages": ["Step status codes: P = pass, C = check, F = fail, - = not run.",
                         "MultiQC (multiqc.html in this folder) gives read, alignment, coverage and Kraken2 plots "
                         "for all samples side by side."],
            "tools": [{"name": "run_summary.py", "version": a.version, "container": "mtb-report:1.1.0"},
                      {"name": "MultiQC", "version": "1.35", "container": "quay.io/biocontainers/multiqc:1.35--pyhdfd78af_2"}],
            "files": {"inputs": [], "outputs": [{"path": "multiqc.html"}, {"path": "run_qc.pdf"}]},
            "results": {"title": "Samples in this run", "columns": ["Sample", "QC", "Category", "Resistance detected",
                                                                    "Indeterminate"],
                        "rows": rows, "widths_mm": [32, 28, 50, 40, 30]},
            "tables": [{"title": "Step status per sample", "columns": ["Sample"] + [x[:2] for x in STEPS[:12]],
                        "rows": status_rows, "widths_mm": [36] + [12] * 12}]}
    json.dump(summ, open("R1_run_qc.summary.json", "w"), indent=2)
    build(summ, yaml.safe_load(open(a.text)), a.out_pdf, [])


if __name__ == "__main__":
    main()
