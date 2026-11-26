#!/usr/bin/env python3
"""Export every variant found for one sample into a variant_export/ folder, with its WHO 2023 classification.

Standalone tool run on finished pipeline results; it does not change any pipeline output.
Usage: export_variants.py SAMPLE_RESULTS_DIR RESOURCES_DIR GENOME_CSQ_TSV
  SAMPLE_RESULTS_DIR  results/<RUN>/<SAMPLE>
  GENOME_CSQ_TSV      genome-wide calls annotated by bcftools csq (written by export_variants.sh)
"""
import csv, os, sys

sdir, res, genome_tsv = sys.argv[1:4]
sample = os.path.basename(sdir.rstrip("/"))
der = os.path.join(res, "who_catalogue", "2023_v2", "derived")
out = os.path.join(sdir, "variant_export")
os.makedirs(out, exist_ok=True)


def tsv(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


here = os.path.dirname(os.path.abspath(__file__))
drugs = tsv(os.path.join(here, "..", "assets", "drugs17.tsv"))
name2code = {}
for d in drugs:
    name2code.setdefault(d["catalogue_drug"], d["code"]) if d["code"] != "PMD" else None
name2code["Para-aminosalicylic acid"] = "PAS"

by_coord, by_name = {}, {}
for r in tsv(os.path.join(der, "who_variants.tsv")):
    if r["drug"] not in name2code:
        continue
    by_coord.setdefault((r["pos"], r["ref"], r["alt"]), []).append(r)
    by_name.setdefault(r["variant"], []).append(r)
gene_rules = {}
for r in tsv(os.path.join(der, "who_gene_rules.tsv")):
    if r["rule"] == "LoF" and r["drug"] in name2code:
        gene_rules.setdefault(r["gene"], []).append(r)
pas = {f"{r['gene']}_{r['mutation']}" for r in tsv(os.path.join(der, "pas_tbdb.tsv"))}
gene_locus = {}  # same rule as the grader: gene LoF rules apply only inside the gene's own locus
for line in open(os.path.join(der, "analysis_regions.bed")):
    f = line.rstrip("\n").split("\t")
    if not f[3].endswith("_wrap"):
        gene_locus.setdefault(f[3], f[4])

evidence = tsv(os.path.join(sdir, "09_resistance", "data", "evidence.tsv"))
used = {}
for e in evidence:
    used.setdefault(e["pos"], set()).add(e["drug"])
mask = []
for line in open(os.path.join(res, "reference", "H37Rv_NC_000962.3", "derived", "low_confidence_mask.bed")):
    f = line.split("\t")
    mask.append((int(f[1]), int(f[2])))


def masked(pos):
    return any(a < pos <= b for a, b in mask)


CALL = {"fixed": "Reported (fixed, AF ≥90%)", "minority": "Reported (minority, AF 10–90%)",
        "below_threshold": "Not reported (below reporting criteria)"}

rows = []
for v in tsv(os.path.join(sdir, "07_annotation", "data", "variants.annotated.tsv")):
    matches = list(by_coord.get((v["pos"], v["ref"], v["alt"]), []))
    for m in by_name.get(v["display_name"], []):
        if m not in matches:
            matches.append(m)
    lof_here = v["lof"] == "yes" and v["csq_locus"] == gene_locus.get(v["gene"])
    lof_rules = gene_rules.get(v["gene"], []) if lof_here else []
    note = (f"gene-disrupting change lies in {v['csq_locus']}, not in {v['gene']}; {v['gene']} LoF rule not "
            "applied" if v["lof"] == "yes" and not lof_here else "")
    who = sorted({f"{name2code[m['drug']]}: {m['grade_label']}" for m in matches} |
                 {f"{name2code[g['drug']]}: {g['gene']}_LoF {g['grade_label']}" for g in lof_rules})
    grades = {int(m["grade"]) for m in matches} | {int(g["grade"]) for g in lof_rules}
    reported_r = sorted(used.get(v["pos"], set()))
    syn = v["csq_consequence"].startswith("synonymous") or (v["snpeff_effect"] or "").startswith("synonymous")
    if reported_r:
        cat = "Resistance - reported"
    elif grades & {1, 2}:
        cat = "Resistance-associated - not reported (" + ("below reporting criteria" if v["call"] == "below_threshold"
                                                          else "see 09_resistance detail") + ")"
    elif v["display_name"] in pas:
        cat = "PAS literature mutation (not WHO-graded)"
    elif 3 in grades:
        cat = "VUS - uncertain significance (WHO grade 3)"
    elif grades & {4, 5}:
        cat = "Not associated with resistance (WHO grade 4/5)"
    elif syn:
        cat = "Synonymous - not in catalogue"
    else:
        cat = "Not in WHO catalogue"
    rows.append({
        "sample": sample, "gene": v["gene"], "locus_tag": v["csq_locus"], "position": v["pos"], "ref": v["ref"],
        "alt": v["alt"], "variant": v["display_name"],
        "who_variant_name": "; ".join(sorted({m["variant"] for m in matches})), "consequence": v["csq_consequence"] or v["snpeff_effect"],
        "protein_change": v["hgvs_p"], "dna_change": v["hgvs_c"], "allele_frequency": v["af"], "depth": v["dp"],
        "alt_reads_fwd": v["alt_fwd"], "alt_reads_rev": v["alt_rev"], "callers": v["callers"],
        "call_status": CALL.get(v["call"], v["call"]), "filter_reason": v["reason"],
        "classification": cat, "drugs_reported_resistant": ", ".join(reported_r),
        "who_2023_grading": "; ".join(who), "loss_of_function": "yes" if lof_here else ("no (other locus)" if v["lof"] == "yes" else "no"),
        "in_low_confidence_mask": "yes" if masked(int(v["pos"])) else "no", "note": note})

ORDER = ["Resistance - reported", "Resistance-associated", "PAS literature", "VUS", "Not associated",
         "Not in WHO catalogue", "Synonymous"]
rows.sort(key=lambda r: (next(i for i, k in enumerate(ORDER) if r["classification"].startswith(k)),
                         r["call_status"].startswith("Not reported"), int(r["position"])))
keys = list(rows[0]) if rows else ["sample"]


def write(name, rs):
    with open(os.path.join(out, f"{sample}.{name}.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, keys, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rs)
    return len(rs)


counts = {"all_resistance_region_variants": write("all_resistance_region_variants", rows),
          "resistance_variants": write("resistance_variants", [r for r in rows if r["classification"].startswith("Resistance")]),
          "vus_variants": write("vus_variants", [r for r in rows if r["classification"].startswith("VUS")]),
          "not_associated_variants": write("not_associated_variants",
                                           [r for r in rows if r["classification"].startswith("Not associated")]),
          "not_in_catalogue_variants": write("not_in_catalogue_variants",
                                             [r for r in rows if r["classification"].startswith(("Not in WHO", "Synonymous", "PAS"))])}

# drug summary (copy of the pipeline's decision, renamed columns)
dr = tsv(os.path.join(sdir, "09_resistance", "data", "drug_results.tsv"))
with open(os.path.join(out, f"{sample}.drug_summary.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["code", "drug", "group", "result", "who_confidence", "evidence_used", "evidence_source",
                "tbprofiler_call", "agrees_with_tbprofiler", "uncertain_variants_for_this_drug"])
    for r in dr:
        w.writerow([r["code"], r["drug"], r["group"], r["call"], r["confidence"], r["detail"], r["evidence_source"],
                    r["tbprofiler"], r["agree"], r["uncertain"]])

# genome-wide calls (bcftools, already csq-annotated by the wrapper)
gw = []
for line in open(genome_tsv):
    pos, ref, alt, qual, dp, dp4, bcsq = line.rstrip("\n").split("\t")
    try:
        a = [int(x) for x in dp4.split(",")]
        af = f"{(a[2] + a[3]) / max(1, sum(a)):.3f}"
    except ValueError:
        af = ""
    gene = cons = change = ""
    if bcsq not in (".", ""):
        f = bcsq.split(",")[0].split("|")
        cons = f[0]
        gene = f[1] if len(f) > 1 else ""
        change = "|".join(x for x in f[5:7] if x) if len(f) > 5 else ""
    gw.append({"sample": sample, "position": pos, "ref": ref, "alt": alt, "gene": gene, "consequence": cons,
               "change": change, "allele_frequency": af, "depth": dp, "quality": qual,
               "in_low_confidence_mask": "yes" if masked(int(pos)) else "no"})
with open(os.path.join(out, f"{sample}.genome_wide_variants.tsv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, list(gw[0]) if gw else ["sample"], delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(gw)
counts["genome_wide_variants"] = len(gw)

cat_n = {}
for r in rows:
    k = r["classification"].split(" (")[0].split(" - ")[0]
    cat_n[k] = cat_n.get(k, 0) + 1
with open(os.path.join(out, "README.txt"), "w") as fh:
    fh.write(f"""Variant export for sample {sample}
Produced by tools/export_variants.py from the pipeline results in this sample folder (no pipeline output changed).
Resistance-gene regions = the 70 WHO 2023 catalogue + PAS gene regions analysed by the pipeline.
Classification uses the WHO 2023 catalogue (grade 1-2 = resistance, 3 = uncertain significance / VUS,
4-5 = not associated with resistance) and, for PAS only, the tbdb literature list.

Files
  {sample}.drug_summary.tsv                    one row per drug: DETECTED / NOT DETECTED / INDETERMINATE (the pipeline's result)
  {sample}.all_resistance_region_variants.tsv  every variant in the resistance-gene regions, classified ({counts['all_resistance_region_variants']} rows)
  {sample}.resistance_variants.tsv             variants graded 1-2 for any drug ({counts['resistance_variants']} rows)
  {sample}.vus_variants.tsv                    WHO grade 3, uncertain significance ({counts['vus_variants']} rows)
  {sample}.not_associated_variants.tsv         WHO grade 4-5, not associated with resistance ({counts['not_associated_variants']} rows)
  {sample}.not_in_catalogue_variants.tsv       not in the WHO catalogue (incl. synonymous, PAS literature) ({counts['not_in_catalogue_variants']} rows)
  {sample}.genome_wide_variants.tsv            every variant called across the whole genome (bcftools), with gene/effect ({counts['genome_wide_variants']} rows)

Counts by classification (resistance-gene regions): {cat_n}

Key columns
  allele_frequency  share of reads carrying the change (1.0 = all bacteria)
  call_status       Reported (fixed AF>=90% / minority AF 10-90%, both minority callers agree, >=5 reads per strand)
                    or Not reported (below the reporting criteria; shown for completeness, never used for a result)
  classification    Resistance - reported / Resistance-associated - not reported / VUS / Not associated /
                    Not in WHO catalogue / Synonymous / PAS literature
  who_2023_grading  per drug: WHO grade label; LoF changes also list the gene-level <gene>_LoF rule
  in_low_confidence_mask  repeats / PE-PPE regions where short-read calls are less reliable
""")
print(sample, counts)
