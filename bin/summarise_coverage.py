#!/usr/bin/env python3
"""Step 04 summary: genome depth/breadth, resistance-gene coverage, gene drop-out screen and the
two-tier key-position counts (the per-drug Indeterminate call itself is made in step 09)."""
import argparse, gzip
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from mtbstep import Step, grade, key_depth, style_axes, SERIES1, SERIES2, INK, INK2, MUTED

GENOME = 4411532


def per_base(path):
    depth = np.zeros(GENOME, dtype=np.int32)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            _, s, e, d = line.split("\t")
            depth[int(s):int(e)] = int(d)
    return depth


def ref_line(ax, y, label):
    ax.axhline(y, color=MUTED, lw=1, ls=":")
    ax.text(1.0, y, f" {label}", color=INK2, fontsize=7.5, va="center", transform=ax.get_yaxis_transform())


def main():
    ap = argparse.ArgumentParser()
    for k in ("summary", "dist", "perbase", "genes", "gene_thresholds", "key", "container", "sample_id",
              "run_id", "thresholds", "version", "out"):
        ap.add_argument("--" + k, required=True)
    a = ap.parse_args()
    s = Step("04_coverage", a.sample_id, a.run_id, a.version, a.thresholds)
    t = s.t

    mean = next(float(l.split()[3]) for l in open(a.summary) if l.startswith("total"))
    dist = {int(l.split()[1]): float(l.split()[2]) for l in open(a.dist) if l.startswith("total")}
    b10 = 100 * dist.get(10, 0)
    s.metric("Mean depth", "Average number of times each DNA letter was read", round(mean, 1), "x",
             f"FAIL < {t['depth_fail_x']}x; CHECK < {t['depth_pass_x']}x",
             grade(mean, t["depth_fail_x"], t["depth_pass_x"]))
    s.metric("Genome at >= 10x", "Share of the TB genome read well enough", round(b10, 2), "%",
             f">= {t['breadth10_pass_pct']}%", grade(b10, t["breadth10_pass_pct"]))
    if mean < t["depth_target_x"]:
        s.msg(f"Mean depth {mean:.0f}x is below the {t['depth_target_x']}x target: resistance present in "
              "10-20% of bacteria may be missed. This is reported per drug, not as a failure.")

    depth = per_base(a.perbase)
    genes = []
    thr = {}
    with gzip.open(a.gene_thresholds, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            ln = int(f[2]) - int(f[1])
            thr[f[3]] = (100 * int(f[4]) / ln, 100 * int(f[5]) / ln)
    with gzip.open(a.genes, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            st, en, g = int(f[1]), int(f[2]), f[3]
            genes.append({"gene": g, "start": st, "mean": float(f[-1]), "min": int(depth[st:en].min()),
                          "b1": thr[g][0], "b10": thr[g][1]})
    drop = [g for g in genes if g["b1"] < t["gene_dropout_pct"]]
    low10 = [g for g in genes if g["b10"] < 100]
    s.metric("Resistance genes fully covered at >= 10x", "Resistance genes read well enough along their whole length",
             f"{len(genes) - len(low10)} of {len(genes)}", threshold="all genes", status="PASS" if not low10 else "WARN")
    s.metric("Possible gene deletions", "Resistance genes with missing pieces",
             ", ".join(g["gene"] for g in drop) or "none",
             threshold=f"each gene >= {t['gene_dropout_pct']}% covered at >= 1x", status="WARN" if drop else "PASS")
    for g in low10:
        edge = " Expected: the gene sits at the start of the linear reference (the genome is circular), " \
               "so reads thin out at position 1." if g["start"] < 300 else ""
        s.msg(f"{g['gene']}: {g['b10']:.1f}% of bases at >= 10x (minimum depth {g['min']}x).{edge}")
    if drop:
        s.msg("Genes flagged as possible deletions are checked by assembly in step 06.")

    lt_ind, lt_min, drugs_ind = 0, 0, set()
    with gzip.open(a.key, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            m = key_depth(depth, int(f[1]), int(f[2]))      # same rule as step 09
            if m < t["key_indeterminate_x"]:
                lt_ind += 1
                drugs_ind |= set(f[4].split(","))
            elif m < t["key_minority_x"]:
                lt_min += 1
    s.metric(f"Key positions < {t['key_indeterminate_x']}x", "Key resistance positions not read well enough",
             lt_ind, threshold="0 (otherwise the drug is Indeterminate in step 09)",
             status="WARN" if lt_ind else "PASS")
    s.metric(f"Key positions {t['key_indeterminate_x']}-{t['key_minority_x']}x",
             "Key positions where low-level resistance cannot be assessed", lt_min)
    if drugs_ind:
        s.msg("Drugs with at least one key position below threshold: " + ", ".join(sorted(drugs_ind)))

    # genome-wide depth in 10 kb windows
    w = 10000
    win = depth[: GENOME // w * w].reshape(-1, w).mean(axis=1)
    fig, ax = plt.subplots(figsize=(9, 3.0))
    ax.plot(np.arange(len(win)) * w / 1e6, win, color=SERIES1, lw=1.2)
    ref_line(ax, t["depth_pass_x"], f"{t['depth_pass_x']}x")
    ref_line(ax, t["key_indeterminate_x"], f"{t['key_indeterminate_x']}x")
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Position in H37Rv genome (Mb)", fontsize=8, color=INK2)
    ax.set_ylabel("Mean depth (10 kb windows)", fontsize=8, color=INK2)
    ax.set_title("Depth along the genome", fontsize=10, color=INK, loc="left")
    style_axes(ax)
    fig.tight_layout()
    fig.savefig("genome_depth.png", dpi=200)

    genes.sort(key=lambda g: g["start"])
    fig, ax = plt.subplots(figsize=(9, 3.6))
    x = np.arange(len(genes))
    ax.scatter(x, [g["mean"] for g in genes], s=22, color=SERIES1, marker="o", label="Mean depth", zorder=3)
    ax.scatter(x, [g["min"] for g in genes], s=26, color=SERIES2, marker="v", label="Lowest depth", zorder=3)
    ref_line(ax, t["key_minority_x"], f"{t['key_minority_x']}x")
    ref_line(ax, t["key_indeterminate_x"], f"{t['key_indeterminate_x']}x")
    ax.set_xticks(x)
    ax.set_xticklabels([g["gene"] for g in genes], rotation=90, fontsize=6.2)
    ax.set_ylim(bottom=0)
    ax.set_ylabel("Depth", fontsize=8, color=INK2)
    ax.set_title("Resistance genes (genome order)", fontsize=10, color=INK, loc="left")
    style_axes(ax)
    ax.legend(loc="upper right", ncol=2, frameon=False, fontsize=7.5, labelcolor=INK2)
    fig.tight_layout()
    fig.savefig("gene_depth.png", dpi=200)

    s.tool("mosdepth", a.container.split(":")[-1].split("--")[0], a.container,
           "mosdepth genome BAM; mosdepth -n --by resistance_regions.bed -T 1,10,50; "
           "mosdepth -n --by graded_positions.bed")
    for p in ("data/genome.per-base.bed.gz", "data/genes.regions.bed.gz", "data/genes.thresholds.bed.gz",
              "data/key.regions.bed.gz"):
        s.output(p)
    s.write(a.out, mean_depth=f"{mean:.0f}", breadth10=f"{b10:.1f}", n_genes=len(genes),
            n_genes_ok=len(genes) - len(low10), key_low=lt_ind,
            dropout_genes=[g["gene"] for g in drop])


if __name__ == "__main__":
    main()
