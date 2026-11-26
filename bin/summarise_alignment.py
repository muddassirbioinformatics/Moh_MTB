#!/usr/bin/env python3
"""Step 03 summary: mapping rate, proper pairs and duplicates from samtools flagstat/markdup JSON."""
import argparse, json
from mtbstep import Step, grade, read_text


def main():
    ap = argparse.ArgumentParser()
    for k in ("flagstat", "markdup", "command_file", "versions_file", "container", "sample_id", "run_id",
              "thresholds", "version", "out"):
        ap.add_argument("--" + k, required=True)
    a = ap.parse_args()
    s = Step("03_alignment", a.sample_id, a.run_id, a.version, a.thresholds)
    t = s.t
    fs = json.load(open(a.flagstat))["QC-passed reads"]
    md = json.load(open(a.markdup))

    primary = fs.get("primary", fs["total"])
    mapped_pct = 100 * fs.get("primary mapped", fs["mapped"]) / primary
    paired_pct = fs.get("properly paired %", 0.0)
    if isinstance(paired_pct, str):
        paired_pct = float(paired_pct.rstrip("%"))
    dup_pct = 100 * fs.get("primary duplicates", fs.get("duplicates", 0)) / primary

    s.metric("Reads mapped to H37Rv", "Share of reads that match the TB reference", round(mapped_pct, 2), "%",
             f"FAIL < {t['mapped_fail_pct']}%; CHECK < {t['mapped_warn_pct']}%",
             grade(mapped_pct, t["mapped_fail_pct"], t["mapped_warn_pct"]))
    s.metric("Properly paired", "Properly paired", round(paired_pct, 2), "%")
    s.metric("Duplicate reads", "Duplicate reads", round(dup_pct, 2), "%")
    s.metric("Primary reads", "Primary reads", primary)
    s.metric("Optical duplicates", "Optical duplicates", md.get("DUPLICATE PAIR OPTICAL", 0) * 2)
    s.msg("Duplicates are marked, not removed; variant callers in step 05 ignore them.")

    vers = read_text(a.versions_file).splitlines()
    s.tool("BWA-MEM", vers[0].split()[-1] if vers else "?", a.container, read_text(a.command_file))
    s.tool("samtools", vers[1].split()[-1] if len(vers) > 1 else "?", a.container)
    for p in (f"data/{a.sample_id}.bam", f"data/{a.sample_id}.bam.bai", "data/flagstat.json",
              "data/markdup.json", "data/samtools_stats.txt"):
        s.output(p)
    s.write(a.out, mapped_pct=f"{mapped_pct:.1f}", dup_pct=f"{dup_pct:.1f}")


if __name__ == "__main__":
    main()
