#!/usr/bin/env python3
"""Step 05 merge: combine bcftools / LoFreq / freebayes (normalised, atomised) calls in the resistance
regions and classify each allele as fixed, minority or below threshold (thresholds.yaml 05_variants).

Writes data/variants.tsv, data/fixed.vcf and data/minority.vcf (sites for annotation in step 07),
the step summary JSON and af_spectrum.png.
"""
import argparse, bisect, gzip
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mtbstep import Step, style_axes, SERIES1, SERIES2, INK, INK2, MUTED

CHROM, GLEN = "NC_000962.3", 4411532


def records(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            info = dict(x.split("=", 1) if "=" in x else (x, True) for x in f[7].split(";"))
            fmt = dict(zip(f[8].split(":"), f[9].split(":"))) if len(f) > 9 else {}
            yield int(f[1]), f[3], f[4], f[5], f[6], info, fmt


def ints(s):
    return [int(x) if x not in (".", "") else 0 for x in str(s).split(",")]


def load(path, caller):
    out = {}
    for pos, ref, alt, qual, filt, info, fmt in records(path):
        if alt in (".", "*"):
            continue
        d = {"dp": int(info.get("DP", 0)), "filter": filt}
        if caller == "bcftools":
            ad = ints(fmt.get("AD", "0,0"))
            tot = sum(ad)
            d["af"] = ad[-1] / tot if tot else 0.0
            d["fwd"], d["rev"] = ints(fmt.get("ADF", "0,0"))[-1], ints(fmt.get("ADR", "0,0"))[-1]
            d["qual"] = float(qual) if qual != "." else 0.0
            d["dp"] = max(d["dp"], tot)
        elif caller == "lofreq":
            d["af"] = float(info.get("AF", 0))
            dp4 = ints(info.get("DP4", "0,0,0,0"))
            d["fwd"], d["rev"] = dp4[2], dp4[3]
        else:
            ro, ao = ints(info.get("RO", 0))[0], ints(info.get("AO", 0))[0]
            d["af"] = ao / (ro + ao) if ro + ao else 0.0
            d["fwd"], d["rev"] = ints(info.get("SAF", 0))[0], ints(info.get("SAR", 0))[0]
        out[(pos, ref, alt)] = d
    return out


class Regions:
    def __init__(self, bed):
        self.iv = sorted((int(f[1]) + 1, int(f[2]), f[3]) for f in (l.split("\t") for l in open(bed)))
        self.starts = [s for s, _, _ in self.iv]

    def gene(self, pos):
        i = bisect.bisect_right(self.starts, pos)
        for s, e, g in reversed(self.iv[max(0, i - 3):i]):
            if s <= pos <= e:
                return g
        return None


def write_vcf(path, rows, sample):
    with open(path, "w") as o:
        o.write("##fileformat=VCFv4.2\n")
        o.write(f"##contig=<ID={CHROM},length={GLEN}>\n")
        o.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency (merged)">\n')
        o.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">\n')
        o.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        o.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n")
        for r in rows:
            o.write(f"{CHROM}\t{r['pos']}\t.\t{r['ref']}\t{r['alt']}\t.\tPASS\tAF={r['af']:.4f};DP={r['dp']}\tGT\t1\n")


def main():
    ap = argparse.ArgumentParser()
    for k in ("bcftools", "lofreq", "freebayes", "regions", "sample_id", "run_id", "thresholds", "version",
              "containers", "out"):
        ap.add_argument("--" + k, required=True)
    a = ap.parse_args()
    s = Step("05_variants", a.sample_id, a.run_id, a.version, a.thresholds)
    t = s.t
    c_b, c_l, c_f = a.containers.split("|")
    min_af, fixed_af = t["min_af_pct"] / 100, t["fixed_af_pct"] / 100
    regions = Regions(a.regions)

    calls = {c: load(p, c) for c, p in (("bcftools", a.bcftools), ("lofreq", a.lofreq),
                                         ("freebayes", a.freebayes))}
    genome = calls["bcftools"]
    n_snp = sum(1 for (p, r, al) in genome if len(r) == len(al))
    n_ind = len(genome) - n_snp

    rows = []
    for key in sorted(set().union(*calls.values())):
        pos, ref, alt = key
        gene = regions.gene(pos)
        if gene is None:
            continue
        b, l, f = (calls[c].get(key) for c in ("bcftools", "lofreq", "freebayes"))
        af = {c: d["af"] for c, d in (("bcftools", b), ("lofreq", l), ("freebayes", f)) if d}
        dp = max(d["dp"] for d in (b, l, f) if d)
        strand = l or f or b
        fwd, rev = strand["fwd"], strand["rev"]
        callers = ",".join(af)
        if dp < t["min_depth"]:
            call, why = "below_threshold", f"depth {dp} < {t['min_depth']}"
        elif (b and b["af"] >= fixed_af and b["qual"] >= 20) or max(af.values()) >= fixed_af:
            call, why = "fixed", ""
        else:
            agree = l is not None and f is not None and l["af"] >= min_af and f["af"] >= min_af
            maf = (l["af"] + f["af"]) / 2 if agree else max(af.values())
            ok_strand = min(fwd, rev) >= t["min_alt_reads_per_strand"]
            lofreq_pass = l is not None and l["filter"] in ("PASS", ".")
            if maf < min_af:
                call, why = "below_threshold", f"AF {maf:.1%} < {min_af:.0%}"
            elif t["minority_caller_agreement"] and not agree:
                call, why = "below_threshold", "LoFreq and freebayes do not agree"
            elif not ok_strand:
                call, why = "below_threshold", f"alt reads per strand {fwd}/{rev} < {t['min_alt_reads_per_strand']}"
            elif not lofreq_pass:
                call, why = "below_threshold", f"LoFreq filter {l['filter']}"
            else:
                call, why = "minority", ""
        if call == "minority":
            vaf = (l["af"] + f["af"]) / 2
        elif l or f:
            vaf = max(x for c, x in af.items() if c != "bcftools") if (l or f) else af["bcftools"]
        else:
            vaf = af["bcftools"]
        rows.append({"pos": pos, "ref": ref, "alt": alt, "gene": gene, "dp": dp, "af": vaf,
                     "af_bcftools": af.get("bcftools"), "af_lofreq": af.get("lofreq"),
                     "af_freebayes": af.get("freebayes"), "alt_fwd": fwd, "alt_rev": rev,
                     "callers": callers, "call": call, "reason": why})

    with open("data/variants.tsv", "w") as o:
        cols = ["pos", "ref", "alt", "gene", "dp", "af", "af_bcftools", "af_lofreq", "af_freebayes",
                "alt_fwd", "alt_rev", "callers", "call", "reason"]
        o.write("\t".join(cols) + "\n")
        for r in rows:
            o.write("\t".join("" if r[c] is None else f"{r[c]:.4f}" if isinstance(r[c], float) else str(r[c])
                              for c in cols) + "\n")
    fixed = [r for r in rows if r["call"] == "fixed"]
    minority = [r for r in rows if r["call"] == "minority"]
    below = [r for r in rows if r["call"] == "below_threshold"]
    write_vcf("data/fixed.vcf", fixed, a.sample_id)
    write_vcf("data/minority.vcf", minority, a.sample_id)

    s.metric("Genome-wide variants (bcftools)", "DNA changes across the whole genome",
             f"{n_snp} SNPs, {n_ind} indels")
    s.metric("Resistance genes: fixed variants", "Changes in resistance genes found in (nearly) all bacteria",
             len(fixed), threshold=f"AF >= {t['fixed_af_pct']}%")
    s.metric("Resistance genes: minority variants",
             "Changes in resistance genes found in only part of the bacteria", len(minority),
             threshold=f"AF {t['min_af_pct']}-{t['fixed_af_pct']}%, >= {t['min_alt_reads_per_strand']} reads/strand, "
                       "LoFreq + freebayes agree")
    s.metric("Resistance genes: below threshold", "Weak signals not used for the result", len(below))
    s.metric("Variant calling completed", "Variant finding completed", "yes", threshold="all three callers ran",
             status="PASS")
    s.msg("Whether a change causes resistance is decided in step 09 (WHO 2023 catalogue), not here.")
    for r in (fixed + minority)[:40]:
        s.msg(f"{r['gene']} {r['pos']} {r['ref']}>{r['alt']}: {r['call']}, AF {r['af']:.1%}, depth {r['dp']}, "
              f"callers {r['callers']}")

    fig, ax = plt.subplots(figsize=(9, 3.0))
    for grp, col, mk, lab in ((fixed, SERIES1, "o", "Fixed"), (minority, SERIES2, "^", "Minority"),
                              (below, MUTED, "x", "Below threshold")):
        if grp:
            ax.scatter([r["pos"] / 1e6 for r in grp], [100 * r["af"] for r in grp], s=22, color=col,
                       marker=mk, label=lab, zorder=3)
    for y in (t["min_af_pct"], t["fixed_af_pct"]):
        ax.axhline(y, color=MUTED, lw=1, ls=":")
        ax.text(1.0, y, f" {y}%", color=INK2, fontsize=7.5, va="center", transform=ax.get_yaxis_transform())
    if not rows:
        ax.text(0.5, 0.5, "No variants in the resistance genes", transform=ax.transAxes, ha="center",
                color=INK2, fontsize=9)
    ax.set_xlim(0, GLEN / 1e6)
    ax.set_ylim(0, 105)
    ax.set_xlabel("Position in H37Rv genome (Mb)", fontsize=8, color=INK2)
    ax.set_ylabel("Allele frequency (%)", fontsize=8, color=INK2)
    ax.set_title("Variants in the resistance genes", fontsize=10, color=INK, loc="left")
    style_axes(ax)
    if rows:
        ax.legend(loc="center right", frameon=False, fontsize=7.5, labelcolor=INK2)
    fig.tight_layout()
    fig.savefig("af_spectrum.png", dpi=200)

    tag = lambda c: c.split(":")[-1].split("--")[0]
    s.tool("bcftools", tag(c_b), c_b, "bcftools mpileup -a AD,ADF,ADR,DP | bcftools call -m --ploidy 1 -v "
           "(genome-wide); bcftools norm -m -any -a (all callers)")
    s.tool("LoFreq", tag(c_l), c_l, "lofreq indelqual --dindel; lofreq call-parallel --call-indels -l resistance_regions.bed")
    s.tool("freebayes", tag(c_f), c_f, "freebayes -p 1 --pooled-continuous -F 0.05 -C 3 -t resistance_regions.bed")
    for p in ("data/variants.tsv", "data/bcftools.norm.vcf.gz", "data/lofreq.norm.vcf.gz",
              "data/freebayes.norm.vcf.gz", "data/fixed.vcf", "data/minority.vcf"):
        s.output(p)
    s.write(a.out, n_fixed=len(fixed), n_minority=len(minority), n_snp=n_snp, n_ind=n_ind)


if __name__ == "__main__":
    main()
