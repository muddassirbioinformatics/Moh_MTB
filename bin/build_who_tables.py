#!/usr/bin/env python3
"""Build the frozen lookup tables the pipeline uses from the WHO 2023 catalogue (run once per catalogue
version; outputs go to resources/who_catalogue/<version>/derived/ and are added to the manifest).

Usage: build_who_tables.py WHO_DIR GFF3 OUT_DIR

Outputs
  who_variants.tsv      drug, gene, variant, pos, ref, alt, grade, grade_label, criteria, comment  (coordinate rows)
  who_gene_rules.tsv    drug, gene, rule (LoF|deletion), grade, grade_label, comment             (gene-level rules)
  resistance_regions.bed  NC_000962.3, start0, end, gene, locus_tag, strand   (gene body + catalogue span)
  graded_positions.bed    NC_000962.3, start0, end, variant, drugs            (grade 1-2 positions; two-tier rule)
"""
import collections, csv, os, sys

CHROM = "NC_000962.3"
who_dir, gff, out = sys.argv[1:4]
os.makedirs(out, exist_ok=True)

# master file: one row per (drug, variant) with the final grade
master = collections.defaultdict(list)
with open(os.path.join(who_dir, "WHO-UCN-TB-2023.6-eng_catalogue_master_file.txt"), newline="") as fh:
    rd = csv.reader(fh, delimiter="\t")
    head = next(rd)
    ix = {k: head.index(k) for k in ("drug", "gene", "variant", "FINAL CONFIDENCE GRADING",
                                     "Additional grading criteria applied", "Comment")}
    for r in rd:
        label = r[ix["FINAL CONFIDENCE GRADING"]]
        master[r[ix["variant"]]].append((r[ix["drug"]], r[ix["gene"]], int(label[0]), label,
                                         r[ix["Additional grading criteria applied"]], r[ix["Comment"]]))

coords = collections.defaultdict(list)
with open(os.path.join(who_dir, "WHO-UCN-TB-2023.7-eng_genomic_coordinates.txt")) as fh:
    next(fh)
    for line in fh:
        v, chrom, pos, ref, alt = line.rstrip("\n").split("\t")
        coords[v].append((int(pos), ref, alt))

n_var = n_rule = 0
missing = []
with open(os.path.join(out, "who_variants.tsv"), "w") as fv, open(os.path.join(out, "who_gene_rules.tsv"), "w") as fr:
    fv.write("drug\tgene\tvariant\tpos\tref\talt\tgrade\tgrade_label\tcriteria\tcomment\n")
    fr.write("drug\tgene\trule\tgrade\tgrade_label\tcomment\n")
    for v, rows in sorted(master.items()):
        for drug, gene, g, label, crit, com in rows:
            if v.endswith("_LoF") or v.endswith("_deletion"):
                fr.write(f"{drug}\t{gene}\t{v.rsplit('_', 1)[1]}\t{g}\t{label}\t{com}\n")
                n_rule += 1
            elif v in coords:
                for pos, ref, alt in coords[v]:
                    fv.write(f"{drug}\t{gene}\t{v}\t{pos}\t{ref}\t{alt}\t{g}\t{label}\t{crit}\t{com}\n")
                    n_var += 1
            else:
                missing.append(v)

# gene bodies from the frozen GFF, keyed by every name WHO might use
genes = {}
for line in open(gff):
    f = line.rstrip("\n").split("\t")
    if len(f) > 8 and f[2] in ("gene", "pseudogene"):
        a = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
        for k in ("gene", "locus_tag", "old_locus_tag", "Name", "gene_synonym"):
            for name in a.get(k, "").split(","):
                if name:
                    genes.setdefault(name, (int(f[3]), int(f[4]), f[6], a.get("locus_tag", "")))

# catalogue positions far from the gene body (dnaA's upstream lies at the end of the circular genome)
# form a separate "<gene>_wrap" interval instead of stretching the gene across the genome
WRAP = 50000
span = collections.defaultdict(lambda: [10**9, 0])
wrap = collections.defaultdict(lambda: [10**9, 0])
for v, rows in master.items():
    gene = rows[0][1]
    gs, ge = genes[gene][:2]
    for pos, ref, alt in coords.get(v, []):
        sp = wrap if (pos < gs - WRAP or pos > ge + WRAP) else span
        sp[gene][0] = min(sp[gene][0], pos)
        sp[gene][1] = max(sp[gene][1], pos + len(ref) - 1)
with open(os.path.join(out, "resistance_regions.bed"), "w") as fb:
    for gene in sorted({r[0][1] for r in master.values()}):
        s, e, strand, lt = genes[gene]
        lo, hi = min(s, span[gene][0]), max(e, span[gene][1])
        fb.write(f"{CHROM}\t{lo - 1}\t{hi}\t{gene}\t{lt}\t{strand}\n")
        if gene in wrap:
            fb.write(f"{CHROM}\t{wrap[gene][0] - 1}\t{wrap[gene][1]}\t{gene}_wrap\t{lt}\t{strand}\n")

key = collections.defaultdict(set)
for v, rows in master.items():
    drugs = {d for d, _, g, _, _, _ in rows if g in (1, 2)}
    if drugs:
        for pos, ref, alt in coords.get(v, []):
            key[(pos, pos + len(ref) - 1, v)] |= drugs
with open(os.path.join(out, "graded_positions.bed"), "w") as fk:
    for (s, e, v), drugs in sorted(key.items()):
        fk.write(f"{CHROM}\t{s - 1}\t{e}\t{v}\t{','.join(sorted(drugs))}\n")

print(f"variant rows {n_var}; gene rules {n_rule}; regions {len(span)}; grade-1/2 intervals {len(key)}; "
      f"graded variants without coordinates {len(missing)}")
if missing:
    print("examples without coordinates:", missing[:10])
