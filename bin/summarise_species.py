#!/usr/bin/env python3
"""Step 02 summary: Mykrobe species call + Sylph profile must agree on MTBC; Kraken2 is a screen only.

Mykrobe's resistance predictions in the same JSON are deliberately ignored (WHO 2023 grading is step 09).
"""
import argparse, csv, json
from mtbstep import Step, grade

MTBC_TAXID = "77643"


def mykrobe_call(path):
    d = json.load(open(path))
    phy = next(iter(d.values()))["phylogenetics"]

    def top(level):
        calls = phy.get(level, {})
        if not calls:
            return None, None
        name = max(calls, key=lambda k: calls[k].get("percent_coverage", 0))
        return name, calls[name]
    out = {lvl: top(lvl) for lvl in ("phylo_group", "sub_complex", "species")}
    lin = phy.get("lineage", {}).get("lineage", [])
    out["lineage"] = (", ".join(lin) if lin else None, None)
    return out, d


def sylph_hits(path, min_ani):
    hits = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ani = float(r["Adjusted_ANI"])
            if ani >= min_ani:
                name = " ".join(r["Contig_name"].split()[1:3])
                hits.append({"name": name, "ani": ani, "abundance": float(r["Taxonomic_abundance"]),
                             "eff_cov": float(r["Eff_cov"]), "genome": r["Genome_file"].split("/")[-1]})
    return sorted(hits, key=lambda h: -h["abundance"])


def kraken(path):
    """Species outside the MTB complex; MTBC members (indented under taxid 77643) are not 'other'."""
    mtbc, unclass, species, mtbc_indent = 0.0, 0.0, [], None
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        pct, rank, taxid, raw = float(f[0]), f[3], f[4], f[5]
        indent = len(raw) - len(raw.lstrip())
        if mtbc_indent is not None and indent <= mtbc_indent:
            mtbc_indent = None
        if taxid == MTBC_TAXID:
            mtbc, mtbc_indent = pct, indent
        elif rank == "U":
            unclass = pct
        elif rank == "S" and mtbc_indent is None:
            species.append((raw.strip(), pct))
    return mtbc, unclass, sorted(species, key=lambda x: -x[1])


def main():
    ap = argparse.ArgumentParser()
    for k in ("mykrobe", "sylph", "kraken2", "sample_id", "run_id", "thresholds", "version", "containers", "out"):
        ap.add_argument("--" + k, required=True)
    a = ap.parse_args()
    s = Step("02_species", a.sample_id, a.run_id, a.version, a.thresholds)
    t = s.t
    c_myk, c_syl, c_k2 = a.containers.split("|")

    myk, raw = mykrobe_call(a.mykrobe)
    group = myk["phylo_group"][0] or "none"
    species = (myk["species"][0] or "not determined").replace("_", " ")
    myk_mtbc = group == "Mycobacterium_tuberculosis_complex"
    s.metric("Mykrobe: MTB complex", "TB family detected (method 1: Mykrobe)",
             "detected" if myk_mtbc else "not detected", threshold="MTBC required",
             status="PASS" if myk_mtbc else "FAIL")
    s.metric("Mykrobe: species", "TB family member (method 1)", species)

    hits = sylph_hits(a.sylph, t["sylph_min_ani"])
    mtb_hits = [h for h in hits if h["name"].startswith("Mycobacterium tuberculosis")]
    syl_mtbc = bool(mtb_hits)
    s.metric("Sylph: MTB complex", "TB family detected (method 2: Sylph)",
             f"ANI {mtb_hits[0]['ani']:.1f}%" if syl_mtbc else "not detected",
             threshold=f"ANI >= {t['sylph_min_ani']}%", status="PASS" if syl_mtbc else "FAIL")
    s.metric("Methods agree", "Both methods agree", "yes" if myk_mtbc == syl_mtbc else "no",
             threshold="Mykrobe = Sylph", status="PASS" if myk_mtbc == syl_mtbc else "WARN")

    others = [h for h in hits if not h["name"].startswith("Mycobacterium tuberculosis")
              and h["abundance"] >= t["other_min_abundance_pct"]]
    s.metric("Other organisms (Sylph)", "Other organisms found in the sample",
             "; ".join(f"{h['name']} ({h['abundance']:.1f}%)" for h in others) or "none",
             threshold=f"none at ANI >= {t['sylph_min_ani']}% and >= {t['other_min_abundance_pct']}%",
             status="WARN" if others else "PASS")

    k_mtbc, k_un, k_sp = kraken(a.kraken2)
    s.metric("Kraken2: reads classified MTBC", "Share of reads that look like TB (screen)", round(k_mtbc, 2), "%",
             f"CHECK < {t['kraken2_mtbc_warn_pct']}%", grade(k_mtbc, warn_below=t["kraken2_mtbc_warn_pct"]))
    s.metric("Kraken2: unclassified", "Unclassified", round(k_un, 2), "%")
    k_other = [(n, p) for n, p in k_sp if p >= t["other_min_abundance_pct"]]
    s.metric("Kraken2: other species >= 1%", "Other species (screen)",
             "; ".join(f"{n} ({p:.1f}%)" for n, p in k_other[:6]) or "none")

    for lvl in ("sub_complex", "lineage"):
        if myk[lvl][0]:
            s.msg(f"Mykrobe {lvl.replace('_', ' ')}: {myk[lvl][0]} (lineage is called in step 08 by TB-Profiler).")
    s.msg("Mykrobe also outputs resistance predictions from its own (WHO 2021-based) panel; these are not used. "
          "Resistance is graded only against the WHO 2023 catalogue in step 09.")
    s.msg("Observed-vs-expected breadth check for co-infecting organisms (plan v1.1) is not yet implemented; "
          "organisms listed here need manual review.")
    for h in hits[:8]:
        s.msg(f"Sylph hit: {h['name']} - ANI {h['ani']:.2f}%, abundance {h['abundance']:.2f}%, "
              f"effective coverage {h['eff_cov']:.1f}x ({h['genome']})")

    v = raw[next(iter(raw))].get("version", {})
    myk_ver = v.get("mykrobe-predictor", "v0.13.0").lstrip("v") + f" (panel {v.get('panel', '?')})"
    s.tool("Mykrobe", myk_ver, c_myk, "mykrobe predict -S tb --panels_dir resources/mykrobe/0.13.0 (panel 202309)")
    s.tool("Sylph", c_syl.split(":")[-1].split("--")[0], c_syl,
           "sylph profile resources/sylph/gtdb_r226_c200_dbv1/gtdb-r226-c200-dbv1.syldb -1 R1 -2 R2")
    s.tool("Kraken2", c_k2.split(":")[-1].split("--")[0], c_k2,
           "kraken2 --db resources/kraken2/standard8_20241228 --paired")
    for p in ("data/mykrobe.json", "data/sylph.tsv", "data/kraken2.report"):
        s.output(p)
    s.write(a.out, species=species, mtbc_word="detected" if (myk_mtbc and syl_mtbc) else "not confirmed",
            others=", ".join(h["name"] for h in others) or "none")


if __name__ == "__main__":
    main()
