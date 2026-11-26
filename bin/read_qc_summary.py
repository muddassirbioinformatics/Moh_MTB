#!/usr/bin/env python3
"""Step 01: turn fastp JSON into the step summary JSON and the per-cycle quality plot.

Usage: read_qc_summary.py --fastp fastp.json --sample_id ID --run_id RUN --thresholds thresholds.yaml
                          --version V --fastp_version_file F --fastp_container IMG --command_file F
                          --out 01_read_qc.summary.json --plot quality.png
"""
import argparse, datetime, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import yaml

# reference palette slots 1 and 2; "after" is also dashed so colour is not the only cue
BEFORE, AFTER = "#2a78d6", "#eb6834"
INK, INK2, MUTED, GRID, AXIS = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"


def plot(fp, out):
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6), sharey=True)
    for ax, rd in zip(axes, ("read1", "read2")):
        for stage, col, ls, lab in (("before", BEFORE, "-", "Before trimming"),
                                    ("after", AFTER, "--", "After trimming")):
            key = f"{rd}_{stage}_filtering"
            if key not in fp:
                continue
            y = fp[key]["quality_curves"]["mean"]
            ax.plot(range(1, len(y) + 1), y, color=col, lw=2, ls=ls, label=lab)
        ax.axhline(30, color=MUTED, lw=1, ls=":")
        ax.text(1.0, 30, " Q30", color=INK2, fontsize=8, va="center", transform=ax.get_yaxis_transform())
        ax.set_title("Read 1" if rd == "read1" else "Read 2", fontsize=10, color=INK, loc="left")
        ax.set_xlabel("Position in read (cycle)", fontsize=8.5, color=INK2)
        ax.grid(axis="y", color=GRID, lw=0.6)
        ax.tick_params(colors=INK2, labelsize=8, length=0)
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.spines["bottom"].set_color(AXIS)
    axes[0].set_ylabel("Mean base quality (Phred)", fontsize=8.5, color=INK2)
    axes[0].set_ylim(bottom=min(20, axes[0].get_ylim()[0]))
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper right", ncol=2, frameon=False, fontsize=8.5, labelcolor=INK2)
    fig.suptitle("Base quality along the read", x=0.01, ha="left", fontsize=11, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(out, dpi=200)


def main():
    ap = argparse.ArgumentParser()
    for k in ("fastp", "sample_id", "run_id", "thresholds", "version", "fastp_version_file", "fastp_container",
              "command_file", "out", "plot"):
        ap.add_argument("--" + k, required=True)
    a = ap.parse_args()
    fp = json.load(open(a.fastp))
    th = yaml.safe_load(open(a.thresholds))
    t = th["01_read_qc"]
    prov = t.get("provisional", False)
    before, after = fp["summary"]["before_filtering"], fp["summary"]["after_filtering"]
    fr = fp["filtering_result"]

    raw_gb = before["total_bases"] / 1e9
    q30 = after["q30_rate"] * 100
    kept = 100 * after["total_reads"] / before["total_reads"]
    metrics = []

    def add(name, plain, value, unit, threshold, status):
        metrics.append({"name": name, "plain": plain, "value": value, "unit": unit, "threshold": threshold,
                        "status": status, "provisional": prov})

    add("Raw yield", "Amount of DNA sequence produced", round(raw_gb, 2), "Gb",
        f"FAIL < {t['raw_yield_fail_gb']} Gb; CHECK < {t['raw_yield_target_gb']} Gb",
        "FAIL" if raw_gb < t["raw_yield_fail_gb"] else "WARN" if raw_gb < t["raw_yield_target_gb"] else "PASS")
    add("Q30 bases after trimming", "Share of high-quality DNA letters", round(q30, 2), "%",
        f">= {t['q30_after_pass_pct']}%", "PASS" if q30 >= t["q30_after_pass_pct"] else "FAIL")
    add("Reads kept after trimming", "Share of reads kept after cleaning", round(kept, 2), "%",
        f"CHECK < {t['reads_kept_warn_pct']}%", "PASS" if kept >= t["reads_kept_warn_pct"] else "WARN")
    info = [("Raw reads", before["total_reads"], ""), ("Reads after trimming", after["total_reads"], ""),
            ("Mean read length R1 after", after.get("read1_mean_length", 0), "bp"),
            ("Q30 bases before trimming", round(before["q30_rate"] * 100, 2), "%"),
            ("Duplication rate", round(fp.get("duplication", {}).get("rate", 0) * 100, 2), "%"),
            ("Insert size peak", fp.get("insert_size", {}).get("peak", 0), "bp"),
            ("GC content after", round(after["gc_content"] * 100, 2), "%")]
    for name, v, u in info:
        metrics.append({"name": name, "plain": name, "value": v, "unit": u, "threshold": "reported only",
                        "status": "INFO"})
    msgs = [f"Reads removed: low quality {fr.get('low_quality_reads', 0):,}; too many N "
            f"{fr.get('too_many_N_reads', 0):,}; too short {fr.get('too_short_reads', 0):,}."]
    if fp.get("adapter_cutting"):
        ac = fp["adapter_cutting"]
        msgs.append(f"Adapter-trimmed reads: {ac.get('adapter_trimmed_reads', 0):,}; "
                    f"adapter-trimmed bases: {ac.get('adapter_trimmed_bases', 0):,}.")

    st = [m["status"] for m in metrics]
    status = "FAIL" if "FAIL" in st else "WARN" if "WARN" in st else "PASS"
    plot(fp, a.plot)
    summ = {
        "schema": "mtb-step-summary/1", "sample_id": a.sample_id, "run_id": a.run_id, "step": "01_read_qc",
        "status": status, "pipeline_version": a.version,
        "generated": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
        "thresholds_version": th.get("version"),
        "values": {"raw_yield_gb": f"{raw_gb:.2f}", "q30_after_pct": f"{q30:.1f}",
                   "reads_kept_pct": f"{kept:.1f}", **t["fastp"]},
        "metrics": metrics, "messages": msgs,
        "tools": [{"name": "fastp", "version": open(a.fastp_version_file).read().split()[-1],
                   "container": a.fastp_container, "command": open(a.command_file).read().strip()},
                  {"name": "read_qc_summary.py", "version": a.version, "container": "mtb-report:1.1.0"}],
        "files": {"inputs": [], "outputs": [{"path": "data/fastp.json"}, {"path": "data/fastp.html"}]},
    }
    json.dump(summ, open(a.out, "w"), indent=2)


if __name__ == "__main__":
    main()
