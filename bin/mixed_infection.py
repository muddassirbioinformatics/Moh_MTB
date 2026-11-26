#!/usr/bin/env python3
"""Step 10: mixed-strain infection.

1. Heterozygous-site test (MixInfect-style, Sobkowiak et al. 2018 principle): SNP sites outside the
   low-confidence mask with depth >= min_depth where the minor allele is 10-90% and has >= 5 reads.
   >= het_sites_mixed such sites -> mixed. The minor-strain share is estimated as the median minor-allele
   frequency of those sites.
2. Lineage mixture: more than one TB-Profiler lineage at >= 10% of marker reads.
3. OmegaTB rule (HaystackAnalytics; >= 5 minor reads per strand, MAPQ > 45) counted for comparison only.
"""
import argparse, gzip, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from mtbstep import Step, style_axes, SERIES1, INK, INK2, MUTED


def load_mask(bed, glen=4411532):
    m = np.zeros(glen + 1, dtype=bool)
    for line in open(bed):
        f = line.split("\t")
        m[int(f[1]) + 1:int(f[2]) + 1] = True       # 1-based positions
    return m


def het_sites(path, t, mask, per_strand=None):
    out = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            ref, alts = f[3], f[4].split(",")
            if mask[int(f[1])]:
                continue
            if len(ref) != 1 or any(len(x) != 1 for x in alts):
                continue
            fmt = dict(zip(f[8].split(":"), f[9].split(":")))
            ad = [int(x) for x in fmt.get("AD", "0").split(",") if x != "."]
            tot = sum(ad)
            if tot < t["min_depth"] or len(ad) < 2:
                continue
            srt = sorted(range(len(ad)), key=lambda i: -ad[i])
            minor_i = srt[1]
            minor = ad[minor_i]
            maf = minor / tot
            if maf * 100 < t["min_minor_af_pct"] or minor < t["min_minor_reads"]:
                continue
            if per_strand is not None:
                adf = [int(x) for x in fmt.get("ADF", "0").split(",")]
                adr = [int(x) for x in fmt.get("ADR", "0").split(",")]
                if min(adf[minor_i], adr[minor_i]) < per_strand:
                    continue
            out.append((int(f[1]), maf))
    return out


def main():
    ap = argparse.ArgumentParser()
    for k in ("mask", "vcf", "omega_vcf", "tbprofiler", "sample_id", "run_id", "thresholds", "version", "container", "out"):
        ap.add_argument("--" + k, required=True)
    a = ap.parse_args()
    s = Step("10_mixed_infection", a.sample_id, a.run_id, a.version, a.thresholds)
    t = s.t

    mask = load_mask(a.mask)
    het = het_sites(a.vcf, t, mask)
    # dense clusters (several "mixed" sites within cluster_bp) are mis-mapping artefacts: count each once
    indep, last = [], -10**9
    for p, m in sorted(het):
        if p - last > t["cluster_bp"]:
            indep.append((p, m))
        last = p
    n = len(indep)
    mixed_het = n >= t["het_sites_mixed"]
    share = float(np.median([m for _, m in indep])) if indep else 0.0
    omega = het_sites(a.omega_vcf, t, mask, per_strand=t["omega_min_minor_reads_per_strand"])

    tbp = json.load(open(a.tbprofiler))
    lins = [(x["lineage"], x.get("fraction", 100.0)) for x in tbp.get("lineage", [])]
    mains = sorted({l for l, f in lins if "." not in l and f >= 10})
    mixed_lin = len(mains) > 1

    s.metric("Heterozygous sites (outside mask)", "Positions where two different DNA letters are present",
             f"{n} independent ({len(het)} before merging sites within {t['cluster_bp']} bp)",
             threshold=f"mixed if >= {t['het_sites_mixed']} independent sites", status="WARN" if mixed_het else "PASS")
    s.metric("Lineage mixture", "Markers from more than one TB family", ", ".join(mains) if mixed_lin else "no",
             threshold="one main lineage", status="WARN" if mixed_lin else "PASS")
    s.metric("Estimated minor strain share", "Estimated share of the second strain",
             f"{share:.0%}" if mixed_het else "-")
    s.metric("OmegaTB-rule sites (comparison)", "Same test using the OmegaTB rule (for comparison)", len(omega),
             threshold=f">= {t['omega_min_minor_reads_per_strand']} minor reads/strand, MAPQ >= {t['omega_min_mapq']}")
    if mixed_het or mixed_lin:
        s.msg("A mixed infection means resistance present in only one strain may appear as a minority variant "
              "(step 05/09); minority resistance calls should be interpreted with this in mind.")
    s.msg("Method: MixInfect-style heterozygous-site count re-implemented in this pipeline (MixInfect and "
          "QuantTB are not available as validated containers). The threshold is provisional until validation "
          "with in-silico mixtures (V4).")
    for p, m in het[:15]:
        s.msg(f"Heterozygous site {p}: minor allele {m:.1%}")

    fig, ax = plt.subplots(figsize=(9, 2.8))
    if het:
        ax.hist([100 * m for _, m in het], bins=np.arange(10, 52, 2), color=SERIES1, edgecolor="white", linewidth=1)
    else:
        ax.text(0.5, 0.5, "No heterozygous sites outside the mask", transform=ax.transAxes, ha="center",
                color=INK2, fontsize=9)
    ax.set_xlim(10, 50)
    ax.set_xlabel("Minor allele frequency (%)", fontsize=8, color=INK2)
    ax.set_ylabel("Sites", fontsize=8, color=INK2)
    ax.set_title("Heterozygous sites outside the low-confidence mask", fontsize=10, color=INK, loc="left")
    style_axes(ax)
    fig.tight_layout()
    fig.savefig("minor_af.png", dpi=200)

    tag = a.container.split(":")[-1].split("--")[0]
    s.tool("bcftools", tag, a.container, "bcftools mpileup -q 20 (and -q 46 for the OmegaTB rule) -a AD,ADF,ADR "
           "| bcftools call -m -v -V indels; low_confidence_mask.bed applied afterwards")
    s.tool("mixed_infection.py", a.version, "mtb-report:1.1.0")
    s.output("data/het_q20.vcf.gz")
    s.output("data/het_omega.vcf.gz")
    verdict = "detected" if (mixed_het or mixed_lin) else "not detected"
    s.write(a.out, verdict=verdict, n_het=n, share=f"{share:.0%}" if mixed_het else "-",
            mixed=bool(mixed_het or mixed_lin))


if __name__ == "__main__":
    main()
