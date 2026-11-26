#!/usr/bin/env python3
"""Step 08 summary: lineage / sub-lineage / spoligotype from TB-Profiler JSON."""
import argparse, json
from mtbstep import Step


def main():
    ap = argparse.ArgumentParser()
    for k in ("tbprofiler", "sample_id", "run_id", "thresholds", "version", "container", "out"):
        ap.add_argument("--" + k, required=True)
    a = ap.parse_args()
    s = Step("08_lineage", a.sample_id, a.run_id, a.version, a.thresholds)
    t = s.t
    d = json.load(open(a.tbprofiler))
    lin = d.get("lineage", [])
    main_l, sub_l = d.get("main_lineage") or "", d.get("sub_lineage") or ""
    fam = next((x.get("family", "") for x in lin if x.get("lineage") == main_l.split(";")[0]), "")
    mains = sorted({x["lineage"] for x in lin if x["lineage"].count(".") == 0 and
                    x.get("fraction", 100) >= t["mixed_lineage_min_pct"]})

    s.metric("Main lineage", "TB family (lineage)", f"{main_l} ({fam})" if main_l else "not determined",
             threshold="determined", status="PASS" if main_l else "WARN")
    s.metric("Sub-lineage", "Sub-family", sub_l or "-")
    s.metric("Single lineage", "Only one TB family detected", "yes" if len(mains) <= 1 else ", ".join(mains),
             threshold=f"one main lineage at >= {t['mixed_lineage_min_pct']}%",
             status="PASS" if len(mains) <= 1 else "WARN")
    for x in lin:
        s.msg(f"{x['lineage']} ({x.get('family', '')}): {x.get('fraction', 0):.1f}% of lineage-marker reads"
              + (f"; RD {x['rd']}" if x.get("rd") not in (None, "None") else ""))
    sp = d.get("spoligotype") or {}
    if sp:
        binary = sp.get("binary", "").replace("■", "1").replace("□", "0")
        s.metric("Spoligotype", "Traditional strain type (from the genome)",
                 f"{sp.get('family', '?')} SIT {sp.get('SIT', '?')} (octal {sp.get('octal', '?')})")
        s.msg(f"Spoligotype binary: {binary}")
    qc = d.get("qc", {})
    if qc:
        s.msg(f"TB-Profiler QC: {qc.get('percent_reads_mapped', '?')}% reads mapped, genome median depth "
              f"{qc.get('genome_median_depth', '?')}x.")
    s.msg("TB-Profiler also predicts resistance from its own database; that is used only as a cross-check "
          "in step 09, never on the report.")
    db = d.get("pipeline", {}).get("db_version", {})
    s.tool("TB-Profiler", d.get("pipeline", {}).get("software_version", "?"), a.container,
           f"tb-profiler profile --no_trim --spoligotype --db tbdb (commit {db.get('commit', '?')}, {db.get('date', '')})")
    s.output("data/tbprofiler.results.json")
    s.output("data/tbprofiler.results.txt")
    s.write(a.out, lineage=main_l or "not determined", family=fam or "unknown", sub_lineage=sub_l or "-",
            spoligo=sp.get("family", "not determined"), sit=sp.get("SIT", "-"))


if __name__ == "__main__":
    main()
