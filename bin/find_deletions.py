#!/usr/bin/env python3
"""Step 06: confirm deletions in resistance genes from contig-to-reference alignments (PAF with cg:Z).

A deletion is "confirmed" when one contig spans it: either a CIGAR D operation >= min_deletion_bp, or two
blocks of the same contig that are adjacent on the contig but separated on the reference. A gene flagged by
coverage drop-out (step 04) without such evidence stays "unconfirmed" (step 09 treats it conservatively).
Without --paf (no drop-out, assembly skipped) the step reports that no assembly was needed.
"""
import argparse, gzip, json, re
import numpy as np
from mtbstep import Step


def regions(bed):
    out = {}
    for line in open(bed):
        f = line.rstrip("\n").split("\t")
        out[f[3]] = (int(f[1]) + 1, int(f[2]))
    return out


def paf_deletions(paf, min_del, max_qgap, max_del, min_mapq, min_block):
    """Return reference intervals (start, end, evidence, contig) deleted in the assembly.
    Only uniquely placed blocks are used, so repeats (IS6110, PE/PPE) cannot fake a join."""
    dels, blocks = [], {}
    for line in open(paf):
        f = line.rstrip("\n").split("\t")
        q, qs, qe, strand, ts, te = f[0], int(f[2]), int(f[3]), f[4], int(f[7]), int(f[8])
        if int(f[11]) < min_mapq or qe - qs < min_block:
            continue
        blocks.setdefault(q, []).append((qs, qe, ts, te, strand))
        cg = next((x[5:] for x in f[12:] if x.startswith("cg:Z:")), "")
        pos = ts
        for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cg):
            n = int(n)
            if op == "D" and min_del <= n <= max_del:
                dels.append((pos + 1, pos + n, "single contig alignment (CIGAR)", q))
            if op in "MD=XN":
                pos += n
    for q, bl in blocks.items():
        bl.sort()
        for (qs1, qe1, ts1, te1, s1), (qs2, qe2, ts2, te2, s2) in zip(bl, bl[1:]):
            if s1 == s2 == "+" and qs2 - qe1 <= max_qgap and min_del <= ts2 - te1 <= max_del:
                dels.append((te1 + 1, ts2, "contig joined across gap", q))
            elif s1 == s2 == "-" and qs2 - qe1 <= max_qgap and min_del <= ts1 - te2 <= max_del:
                dels.append((te2 + 1, ts1, "contig joined across gap", q))
    return dels


def per_base(path, glen=4411532):
    depth = np.zeros(glen, dtype=np.int32)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            _, s, e, d = line.split("\t")
            depth[int(s):int(e)] = int(d)
    return depth


def n50(path):
    lens, cur = [], 0
    for line in open(path):
        if line.startswith(">"):
            if cur:
                lens.append(cur)
            cur = 0
        else:
            cur += len(line.strip())
    lens.append(cur)
    lens.sort(reverse=True)
    acc, half = 0, sum(lens) / 2
    for L in lens:
        acc += L
        if acc >= half:
            return len(lens), sum(lens), L
    return len(lens), sum(lens), 0


def main():
    ap = argparse.ArgumentParser()
    for k in ("coverage_summary", "regions", "sample_id", "run_id", "thresholds", "version", "containers", "out"):
        ap.add_argument("--" + k, required=True)
    ap.add_argument("--perbase", required=True)
    ap.add_argument("--paf")
    ap.add_argument("--contigs")
    a = ap.parse_args()
    s = Step("06_deletions", a.sample_id, a.run_id, a.version, a.thresholds)
    t = s.t
    c_sh, c_mm = a.containers.split("|")
    cov = json.load(open(a.coverage_summary))["values"]
    flagged = cov.get("dropout_genes", [])
    genome_mean = float(cov["mean_depth"])
    depth = per_base(a.perbase)
    reg = regions(a.regions)
    rows = []

    if a.paf:
        n_ctg, total, ctg_n50 = n50(a.contigs)
        s.msg(f"Assembly: {n_ctg} contigs, {total:,} bp, N50 {ctg_n50:,} bp.")
        dels = paf_deletions(a.paf, t["min_deletion_bp"], t["max_query_gap_bp"], t["max_deletion_bp"],
                             t["min_mapq"], t["min_block_bp"])
        # a real deletion has (almost) no reads; assembly-only gaps with normal depth are artefacts
        kept = []
        for ds, de, ev, ctg in dels:
            d = float(depth[ds - 1:de].mean())
            if d <= t["max_depth_fraction"] * genome_mean:
                kept.append((ds, de, f"{ev}; reads {d:.0f}x", ctg))
            elif any(gs <= de and ds <= ge for gs, ge in reg.values()):
                s.msg(f"Assembly gap {ds}-{de} ({de - ds + 1} bp) not reported: reads cover it at {d:.0f}x "
                      f"(genome {genome_mean:.0f}x), so it is an assembly artefact.")
        dels = kept
        for gene, (gs, ge) in reg.items():
            for ds, de, ev, ctg in dels:
                ov = min(ge, de) - max(gs, ds) + 1
                if ov > 0:
                    rows.append({"gene": gene, "start": ds, "end": de, "length": de - ds + 1,
                                 "pct_gene_deleted": round(100 * ov / (ge - gs + 1), 1), "evidence": f"{ev}; {ctg}",
                                 "status": "confirmed"})
        for g in flagged:
            if not any(r["gene"] == g for r in rows):
                rows.append({"gene": g, "start": "", "end": "", "length": "", "pct_gene_deleted": "",
                             "evidence": "coverage drop-out only; assembly does not resolve the region",
                             "status": "unconfirmed"})
        s.tool("shovill (SPAdes)", c_sh.split(":")[-1].split("--")[0], c_sh,
               "shovill --gsize 4.4M --depth 100 (reads subsampled to ~100x)")
        s.tool("minimap2", c_mm.split(":")[-1].split("--")[0], c_mm, "minimap2 -cx asm5 H37Rv.fasta contigs.fa")
        s.output("data/contigs.fa")
        s.output("data/contigs_vs_H37Rv.paf")
    else:
        s.msg("No resistance gene showed coverage drop-out in step 04, so de novo assembly was not needed. "
              "Small insertions/deletions are detected by the variant callers in step 05.")

    with open("data/deletions.tsv", "w") as o:
        cols = ["gene", "start", "end", "length", "pct_gene_deleted", "evidence", "status"]
        o.write("\t".join(cols) + "\n")
        for r in rows:
            o.write("\t".join(str(r[c]) for c in cols) + "\n")
    s.output("data/deletions.tsv")

    conf = [r for r in rows if r["status"] == "confirmed"]
    unconf = [r for r in rows if r["status"] == "unconfirmed"]
    s.metric("Genes flagged by coverage (step 04)", "Genes that looked like they might be missing",
             ", ".join(flagged) or "none")
    s.metric("Assembly run", "Extra check by rebuilding the genome", "yes" if a.paf else "not needed")
    s.metric("Confirmed deletions in resistance genes", "Missing pieces of resistance genes (confirmed)",
             "; ".join(f"{r['gene']} {r['start']}-{r['end']} ({r['length']} bp)" for r in conf) or "none",
             threshold="reported finding", status="WARN" if conf else "PASS")
    s.metric("Unconfirmed drop-outs", "Possible missing pieces that could not be confirmed",
             ", ".join(r["gene"] for r in unconf) or "none", threshold="none",
             status="WARN" if unconf else "PASS")
    s.write(a.out, n_conf=len(conf), n_unconf=len(unconf), assembled="yes" if a.paf else "no",
            deleted=", ".join(sorted({r["gene"] for r in conf})) or "none")


if __name__ == "__main__":
    main()
