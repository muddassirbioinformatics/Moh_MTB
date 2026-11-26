#!/usr/bin/env python3
"""Step 07: add bcftools csq consequences (used for loss-of-function rules) and SnpEff HGVS names
(display only) to the step 05 variant table -> data/variants.annotated.tsv."""
import argparse, collections, csv
from mtbstep import Step

LOF = ("frameshift", "stop_gained", "start_lost")
RNA_GENES = ("rrs", "rrl")


def info_field(line, key):
    for x in line.split("\t")[7].split(";"):
        if x.startswith(key + "="):
            return x[len(key) + 1:]
    return None


def load_csq(paths):
    out = {}
    for p in paths:
        recs = [l.rstrip("\n") for l in open(p) if not l.startswith("#")]
        by_pos = {}
        for l in recs:
            f = l.split("\t")
            by_pos[f[1]] = info_field(l, "BCSQ")
        for l in recs:
            f = l.split("\t")
            b = by_pos[f[1]]
            if b and b.startswith("@"):          # joint codon call reported at another record
                b = by_pos.get(b[1:].split(",")[0], b)
            out[(int(f[1]), f[3], f[4])] = b or ""
    return out


def load_snpeff(paths):
    out = {}
    for p in paths:
        for l in open(p):
            if l.startswith("#"):
                continue
            f = l.split("\t")
            out[(int(f[1]), f[3], f[4])] = [a.split("|") for a in (info_field(l, "ANN") or "").split(",") if a]
    return out


def pick_ann(anns, gene):
    for a in anns:
        if len(a) > 10 and gene in (a[3], a[4]):
            return a
    for a in anns:
        if len(a) > 10 and "intergenic" not in a[1]:
            return a
    return anns[0] if anns else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variants", required=True)
    ap.add_argument("--csq", nargs=2, required=True)
    ap.add_argument("--snpeff", nargs=2, required=True)
    ap.add_argument("--regions", required=True)
    for k in ("sample_id", "run_id", "thresholds", "version", "containers", "out"):
        ap.add_argument("--" + k, required=True)
    a = ap.parse_args()
    s = Step("07_annotation", a.sample_id, a.run_id, a.version, a.thresholds)
    c_b, c_s = a.containers.split("|")
    csq, snp = load_csq(a.csq), load_snpeff(a.snpeff)

    # analysis windows overlap for neighbouring genes (ethA/ethR, embA/embB/embC, Rv0678/mmpS5/mmpL5,
    # fbiA/fbiB, lpqB/mtrB): a change inside a gene is assigned to the gene csq places it in
    locus_gene = {}
    for line in open(a.regions):
        f = line.rstrip("\n").split("\t")
        if not f[3].endswith("_wrap"):
            locus_gene.setdefault(f[4], f[3])
    rows = list(csv.DictReader(open(a.variants), delimiter="\t"))
    used = [r for r in rows if r["call"] in ("fixed", "minority")]
    kinds = collections.Counter()
    missing = 0
    for r in rows:
        key = (int(r["pos"]), r["ref"], r["alt"])
        b = csq.get(key, "")
        parts = b.split(",")[0].split("|") if b else []
        r["csq_consequence"] = parts[0] if parts else ""
        r["csq_change"] = parts[5] if len(parts) > 5 else ""
        r["csq_locus"] = parts[2] if len(parts) > 2 else ""
        if r["csq_locus"] in locus_gene and locus_gene[r["csq_locus"]] != r["gene"]:
            r["window_gene"], r["gene"] = r["gene"], locus_gene[r["csq_locus"]]
        ann = pick_ann(snp.get(key, []), r["gene"])
        r["snpeff_effect"] = ann[1] if ann else ""
        r["hgvs_c"] = ann[9] if ann else ""
        r["hgvs_p"] = ann[10] if ann else ""
        if r["gene"] in RNA_GENES:
            name = f"n.{r['hgvs_c'][2:]}" if r["hgvs_c"] else ""
        else:
            name = r["hgvs_p"] or r["hgvs_c"]
        r["display_name"] = f"{r['gene']}_{name}" if name else f"{r['gene']} {r['pos']} {r['ref']}>{r['alt']}"
        r["lof"] = "yes" if any(k in r["csq_consequence"] for k in LOF) else "no"
        if r["call"] in ("fixed", "minority"):
            if not (r["csq_consequence"] or r["snpeff_effect"]):
                missing += 1
            kinds[(r["csq_consequence"] or r["snpeff_effect"] or "unannotated").split("&")[0]] += 1

    cols = list(rows[0].keys()) if rows else ["pos", "ref", "alt", "gene", "call"]
    with open("data/variants.annotated.tsv", "w") as o:
        w = csv.DictWriter(o, fieldnames=cols + [c for c in ("window_gene", "csq_consequence", "csq_change", "csq_locus", "snpeff_effect",
                                                              "hgvs_c", "hgvs_p", "display_name", "lof")
                                                if c not in cols], delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    s.metric("Variants annotated", "Changes described in terms of their effect on the gene",
             f"{len(used) - missing} of {len(used)}", threshold="all fixed/minority variants",
             status="PASS" if not missing else "WARN")
    s.metric("Gene-disrupting changes (LoF)", "Changes that break a gene",
             sum(1 for r in used if r["lof"] == "yes"))
    for k, n in kinds.most_common():
        s.metric(f"Consequence: {k}", f"Type: {k.replace('_', ' ')}", n)
    for r in used[:40]:
        s.msg(f"{r['display_name']} ({r['call']}, AF {float(r['af']):.1%}): {r['csq_consequence'] or r['snpeff_effect']}"
              + (" - gene-disrupting" if r["lof"] == "yes" else ""))
    s.msg("Consequences from bcftools csq decide loss-of-function; SnpEff names are for display. "
          "rRNA (rrs, rrl) changes are named n.<position> and matched by coordinate in step 09.")
    tag = lambda c: c.split(":")[-1].split("--")[0]
    s.tool("bcftools csq", tag(c_b), c_b, "bcftools csq -C 11 -p a (fixed); -l (minority) with derived/H37Rv.csq.gff3")
    s.tool("SnpEff", tag(c_s), c_s, "snpEff -ud 200 H37Rv (local database built from the frozen GenBank)")
    for p in ("data/variants.annotated.tsv", "data/fixed.csq.vcf", "data/minority.csq.vcf",
              "data/fixed.snpeff.vcf", "data/minority.snpeff.vcf"):
        s.output(p)
    s.write(a.out, n_used=len(used), n_lof=sum(1 for r in used if r["lof"] == "yes"))


if __name__ == "__main__":
    main()
