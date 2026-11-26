#!/usr/bin/env python3
"""Step 09: grade resistance per drug against the WHO 2023 catalogue (exact coordinate match, as WHO
instructs), plus the catalogue's gene-level LoF/deletion rules, codon-level matching for multi-base
changes, the two epistasis rules, pretomanid via the WHO DLM-PMD cross-resistance rule and PAS via the
tbdb literature list (labelled lower evidence). TB-Profiler's own calls are compared as a cross-check.

Outputs data/drug_results.tsv, data/evidence.tsv and the step summary JSON.
"""
import argparse, collections, csv, gzip, json, re
import numpy as np
from mtbstep import Step, key_depth

AA3 = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu", "G": "Gly",
       "H": "His", "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser",
       "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val", "*": "*"}
GLEN = 4411532
ABROGATE = re.compile(r"Abrogates effect of genetically linked AwR and AwRI (\S+) mutations")


def tsv(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def per_base(path):
    depth = np.zeros(GLEN, dtype=np.int32)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            _, s, e, d = line.split("\t")
            depth[int(s):int(e)] = int(d)
    return depth


def main():
    ap = argparse.ArgumentParser()
    for k in ("variants", "deletions", "coverage_summary", "perbase", "derived", "drugs", "tbprofiler",
              "sample_id", "run_id", "thresholds", "version", "out"):
        ap.add_argument("--" + k, required=True)
    a = ap.parse_args()
    s = Step("09_resistance", a.sample_id, a.run_id, a.version, a.thresholds)
    t = s.t
    D = a.derived

    # ---- reference tables
    by_coord = collections.defaultdict(list)
    by_name = collections.defaultdict(list)
    drug_genes = collections.defaultdict(set)
    for r in tsv(f"{D}/who_variants.tsv"):
        r["grade"] = int(r["grade"])
        by_coord[(int(r["pos"]), r["ref"], r["alt"])].append(r)
        by_name[r["variant"]].append(r)
        drug_genes[r["drug"]].add(r["gene"])
    rules = collections.defaultdict(list)
    for r in tsv(f"{D}/who_gene_rules.tsv"):
        r["grade"] = int(r["grade"])
        rules[(r["gene"], r["rule"])].append(r)
        drug_genes[r["drug"]].add(r["gene"])
    pas = tsv(f"{D}/pas_tbdb.tsv")
    for r in pas:
        drug_genes[r["drug"]].add(r["gene"])
    drugs = tsv(a.drugs)
    gene_locus = {}
    for line in open(f"{D}/analysis_regions.bed"):
        f = line.rstrip("\n").split("\t")
        gene_locus.setdefault(f[3], f[4])

    variants = [r for r in tsv(a.variants) if r["call"] in ("fixed", "minority")]
    deletions = [r for r in tsv(a.deletions)]
    cov = json.load(open(a.coverage_summary))["values"]
    unconfirmed = set(cov.get("dropout_genes", [])) - {r["gene"] for r in deletions if r["status"] == "confirmed"}

    # ---- evidence per catalogue drug: dicts with variant, gene, grade, label, af, call, source, comment
    ev = collections.defaultdict(list)

    def add(drug, row, v, source):
        name = v["display_name"]
        if source.startswith("WHO 2023 (") and row.get("variant"):
            name = row["variant"]
            if row.get("comment", "").startswith("Alias "):
                name += f" (alias {row['comment'][6:].split('.', 1)[0]}.{row['comment'][6:].split('.', 1)[1].split('.')[0]})"
        ev[drug].append({"variant": name, "who_variant": row.get("variant", ""), "gene": row["gene"],
                         "grade": row["grade"], "label": row["grade_label"], "af": float(v["af"]),
                         "call": v["call"], "source": source, "comment": row.get("comment", ""), "pos": v["pos"]})

    # codon groups: >= 2 fixed changes csq reports as one missense codon change -> match by name instead
    codon = collections.defaultdict(list)
    for v in variants:
        m = re.match(r"^(\d+)([A-Z*])>\1([A-Z*])$", v.get("csq_change", ""))
        if v["call"] == "fixed" and m and v["csq_consequence"].startswith("missense"):
            codon[(v["gene"], v["csq_change"])].append(v)
    codon_rows = {id(x) for grp in codon.values() if len(grp) > 1 for x in grp}
    for (gene, change), grp in codon.items():
        if len(grp) > 1:
            m = re.match(r"^(\d+)([A-Z*])>\1([A-Z*])$", change)
            name = f"{gene}_p.{AA3[m.group(2)]}{m.group(1)}{AA3[m.group(3)]}"
            v = dict(grp[0], display_name=name)
            for row in by_name.get(name, []):
                add(row["drug"], row, v, "WHO 2023 (codon-level match)")
            if name not in by_name:
                v["unmatched"] = True

    for v in variants:
        if id(v) in codon_rows:
            continue
        key = (int(v["pos"]), v["ref"], v["alt"])
        hits = by_coord.get(key, [])
        seen = set()
        for row in hits:
            if (row["drug"], row["variant"]) not in seen:
                seen.add((row["drug"], row["variant"]))
                add(row["drug"], row, v, "WHO 2023 (coordinate match)")
        lof_here = v["lof"] == "yes" and v.get("csq_locus") == gene_locus.get(v["gene"])
        if v["lof"] == "yes" and not lof_here:
            s.msg(f"{v['display_name']}: gene-disrupting in {v.get('csq_locus')}, not in {v['gene']} "
                  f"({gene_locus.get(v['gene'])}); {v['gene']} LoF rule not applied.")
        if lof_here:
            for row in rules.get((v["gene"], "LoF"), []):
                add(row["drug"], dict(row, variant=f"{row['gene']}_LoF"), v, "WHO 2023 gene rule (LoF)")
        v["unmatched"] = not hits and not lof_here
        for p in pas:
            if p["gene"] == v["gene"] and p["mutation"] in (v.get("hgvs_p"), v.get("hgvs_c")):
                add(p["drug"], {"gene": p["gene"], "grade": 2, "grade_label": "Literature (tbdb, not WHO-graded)",
                                "variant": f"{p['gene']}_{p['mutation']}", "comment": p["comment"]}, v,
                    "tbdb literature")

    for d in (r for r in deletions if r["status"] == "confirmed"):
        v = {"display_name": f"{d['gene']} deletion {d['start']}-{d['end']} ({d['length']} bp)", "af": 1.0,
             "call": "fixed", "pos": d["start"]}
        kinds = ["LoF"] + (["deletion"] if float(d["pct_gene_deleted"]) >= t["delete_rule_min_pct"] else [])
        for k in kinds:
            for row in rules.get((d["gene"], k), []):
                add(row["drug"], dict(row, variant=f"{row['gene']}_{k}"), v, f"WHO 2023 gene rule ({k}, assembly-confirmed)")
        if float(d["pct_gene_deleted"]) >= t["delete_rule_min_pct"]:
            for p in pas:
                if p["gene"] == d["gene"] and p["mutation"] == "transcript_ablation":
                    add(p["drug"], {"gene": p["gene"], "grade": 2, "grade_label": "Literature (tbdb, not WHO-graded)",
                                    "variant": f"{p['gene']}_deletion", "comment": p["comment"]}, v, "tbdb literature")

    # epistasis: a fixed abrogating change cancels fixed AwR/AwRI changes in the named gene
    abrogated_genes = {}
    for drug, es in ev.items():
        for e in es:
            m = ABROGATE.search(e["comment"])
            if m and e["call"] == "fixed":
                abrogated_genes.setdefault((drug, m.group(1)), e["variant"])

    # ---- coverage at graded positions per drug
    depth = per_base(a.perbase)
    low, mid = collections.defaultdict(list), collections.defaultdict(list)
    for line in open(f"{D}/analysis_key_positions.bed"):
        f = line.rstrip("\n").split("\t")
        dmin = key_depth(depth, int(f[1]), int(f[2]))
        for drug in f[4].split(","):
            if dmin < t["indeterminate_depth_x"]:
                low[drug].append(f"{f[3]} ({dmin}x)")
            elif dmin < t["minority_assessable_depth_x"]:
                mid[drug].append(f[3])

    # ---- TB-Profiler cross-check
    tbp = json.load(open(a.tbprofiler))
    tbp_r = set()
    for v in tbp.get("dr_variants", []):
        for dd in v.get("drugs", []):
            tbp_r.add(dd.get("drug"))

    results, ev_rows, ev_full, uncertain_rows = [], [], [], []
    n_ind = n_det = n_disc = 0
    for dr in drugs:
        cat = dr["catalogue_drug"]
        es = ev.get(cat, [])
        if dr["code"] == "PMD":
            es = [e for e in es if "DLM-PMD" in e["comment"]]
        live, abro = [], []
        for e in es:
            gene_abro = abrogated_genes.get((cat, e["gene"]))
            if e["grade"] in (1, 2) and gene_abro and e["call"] == "fixed":
                abro.append((e, gene_abro))
            elif e["grade"] in (1, 2):
                live.append(e)
        # novel / uncertain changes in this drug's genes
        unc = [e for e in es if e["grade"] == 3]
        genes = drug_genes.get(cat, set()) if dr["code"] != "PMD" else {e["gene"] for e in es}
        for v in variants:
            if v.get("unmatched") and v["gene"] in genes and "synonymous" not in (v.get("csq_consequence", "") +
                                                                                 v.get("snpeff_effect", "")):
                unc.append({"variant": v["display_name"], "grade": None, "label": "Not in catalogue",
                            "af": float(v["af"]), "call": v["call"]})
        drug_low = low.get(cat, []) if dr["code"] != "PMD" else []
        gene_unconf = [g for g in unconfirmed if g in genes]

        if live:
            best = min(e["grade"] for e in live)
            minority = all(e["call"] == "minority" for e in live)
            call = "DETECTED"
            conf = "Assoc w R" if best == 1 else "Assoc w R - Interim"
            if dr["evidence"] == "tbdb_literature":
                conf = "Literature only (not WHO-graded)"
            detail = "; ".join(sorted({f"{e['variant']} ({e['af']:.0%})" for e in live}))
            if minority:
                detail += " - minority population"
            n_det += 1
        elif drug_low or gene_unconf:
            call, conf = "INDETERMINATE", ""
            detail = ("key positions below "
                      f"{t['indeterminate_depth_x']}x: {', '.join(drug_low[:5])}" if drug_low else "") + \
                     (f"; unconfirmed gene drop-out: {', '.join(gene_unconf)}" if gene_unconf else "")
            n_ind += 1
        else:
            call, conf, detail = "NOT DETECTED", "", "No resistance-associated change found"
        if call != "DETECTED" and mid.get(cat) and dr["code"] != "PMD":
            detail += f"; low-level resistance not assessable at {len(mid[cat])} position(s) (10-50x)"
        if dr["evidence"] == "tbdb_literature":
            detail += "; PAS is not in the WHO 2023 catalogue - literature list only"
        if dr["code"] == "PMD":
            detail += "; inferred from WHO DLM-PMD cross-resistance (LoF in ddn/fbiA/fbiB/fbiC/fgd1/Rv2983) only"
        for e, by in abro:
            detail += f"; {e['variant']} cancelled by {by} (WHO epistasis rule)"
        tb = "R" if dr["tbprofiler_name"] in tbp_r else "S"
        ours = "R" if call == "DETECTED" else ("?" if call == "INDETERMINATE" else "S")
        agree = "yes" if ours == tb else ("n/a" if ours == "?" else "no")
        n_disc += agree == "no"
        results.append({"code": dr["code"], "drug": dr["name"], "group": dr["group"], "call": call,
                        "confidence": conf, "detail": detail.strip("; "), "evidence_source": dr["evidence"],
                        "tbprofiler": tb, "agree": agree,
                        "uncertain": "; ".join(f"{u['variant']} ({u['af']:.0%}, {u['label']})" for u in unc)})
        for e in live + [x for x, _ in abro]:
            ev_rows.append([dr["code"], e["variant"], e["who_variant"], e["label"], f"{e['af']:.1%}", e["call"],
                            e["source"]])
            ev_full.append([dr["code"], e["variant"], e["who_variant"], e["label"], f"{e['af']:.4f}", e["call"],
                            e["source"], str(e["pos"]), "abrogated" if any(e is x for x, _ in abro) else "used"])
        for u in unc:
            uncertain_rows.append([dr["code"], u["variant"], u["label"], f"{u['af']:.1%}", u["call"]])

    with open("data/drug_results.tsv", "w") as o:
        w = csv.DictWriter(o, fieldnames=list(results[0]), delimiter="\t")
        w.writeheader()
        w.writerows(results)
    with open("data/evidence.tsv", "w") as o:
        o.write("drug\tvariant\twho_variant\tgrade\taf\tcall\tsource\tpos\tstatus\n")
        o.writelines("\t".join(r) + "\n" for r in ev_full)

    s.d["results"] = {"title": "Result per drug", "columns": ["Drug", "Resistance", "Details"],
                      "badge_col": 1, "widths_mm": [38, 34, 108],
                      "rows": [[f"{r['drug']} ({r['code']})", r["call"],
                                (r["detail"] if r["call"] != "NOT DETECTED" else "") or "-"] for r in results]}
    s.d["tables"] = [
        {"title": "Resistance evidence (grade 1-2 matches)", "columns": ["Drug", "Variant", "WHO variant", "Grade",
                                                                       "AF", "Call", "Source"],
         "rows": ev_rows, "widths_mm": [12, 34, 34, 30, 12, 16, 42]},
        {"title": "Changes of uncertain significance / not in catalogue (do not change the result)",
         "columns": ["Drug", "Variant", "Category", "AF", "Call"], "rows": uncertain_rows,
         "widths_mm": [14, 60, 56, 20, 30]},
        {"title": "Cross-check against TB-Profiler (not used for the result)",
         "columns": ["Drug", "This pipeline", "TB-Profiler", "Agree"],
         "rows": [[r["code"], r["call"], "resistant" if r["tbprofiler"] == "R" else "no call", r["agree"]]
                  for r in results], "widths_mm": [20, 60, 60, 40]},
    ]
    s.metric("Drugs assessed", "Drugs with a clear result", f"{len(results) - n_ind} of {len(results)}",
             threshold="all (otherwise Indeterminate)", status="PASS" if not n_ind else "WARN")
    s.metric("Agreement with TB-Profiler", "Second-opinion software agrees",
             "all drugs" if not n_disc else f"{n_disc} drug(s) differ",
             threshold="no discordance (review otherwise)", status="PASS" if not n_disc else "WARN")
    s.metric("Epistasis rules applied", "Cancelling mutations found",
             ", ".join(f"{g} ({d})" for d, g in abrogated_genes) or "none")
    s.msg("Grading: WHO 2023 catalogue (2nd edition), exact (position, ref, alt) match after bcftools "
          "normalisation, as instructed by WHO; grade 1 'Assoc w R' and grade 2 'Assoc w R - Interim' count "
          "as resistance. Grade 3 and changes not in the catalogue are listed but do not change the result.")
    s.msg("Pretomanid is inferred only from the WHO DLM-PMD cross-resistance rule; PAS uses the tbdb "
          "literature list (lower evidence, not WHO-graded).")
    s.tool("grade_resistance.py", a.version, "mtb-report:1.1.0",
           "WHO 2023 derived tables: who_variants.tsv, who_gene_rules.tsv, analysis_key_positions.bed, pas_tbdb.tsv")
    for p in ("data/drug_results.tsv", "data/evidence.tsv"):
        s.output(p)
    det = [f"{r['drug']}" for r in results if r["call"] == "DETECTED"]
    ind = [f"{r['drug']}" for r in results if r["call"] == "INDETERMINATE"]
    s.write(a.out, detected=", ".join(det) or "none", n_detected=len(det), indeterminate=", ".join(ind) or "none",
            n_not=len(results) - len(det) - len(ind), n_drugs=len(results))


if __name__ == "__main__":
    main()
