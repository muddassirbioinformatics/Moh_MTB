#!/usr/bin/env python3
"""Combine the WHO 2023 targets with the secondary (non-WHO) evidence used for PAS into the files the
pipeline actually analyses. WHO tables stay untouched; PAS evidence is labelled as literature (tbdb).

Usage: build_analysis_targets.py WHO_DERIVED_DIR TBDB_DIR GFF3
Writes into WHO_DERIVED_DIR: analysis_regions.bed, analysis_key_positions.bed, pas_tbdb.tsv
"""
import json, os, sys

CHROM = "NC_000962.3"
PAS = "Para-aminosalicylic acid"
who, tbdb, gff = sys.argv[1:4]

genes = {}
for line in open(gff):
    f = line.rstrip("\n").split("\t")
    if len(f) > 8 and f[2] == "gene":
        a = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
        genes[a.get("locus_tag")] = (int(f[3]), int(f[4]), f[6], a.get("gene", a.get("locus_tag")))

muts = json.load(open(os.path.join(tbdb, "mutations.json")))
rows, pas_genes, key = [], set(), set()
for lt, ms in muts.items():
    for m, v in ms.items():
        for an in v["annotations"]:
            if an.get("drug") == "para-aminosalicylic_acid":
                name = genes[lt][3]
                rows.append((PAS, name, lt, m, ",".join(map(str, v.get("genome_positions") or [])),
                             an.get("comment", "")))
                pas_genes.add(lt)
                for p in v.get("genome_positions") or []:
                    key.add((p, name))

with open(os.path.join(who, "pas_tbdb.tsv"), "w") as o:
    o.write("drug\tgene\tlocus_tag\tmutation\tpositions\tcomment\n")
    for r in sorted(rows):
        o.write("\t".join(r) + "\n")

with open(os.path.join(who, "analysis_regions.bed"), "w") as o:
    o.writelines(open(os.path.join(who, "resistance_regions.bed")))
    for lt in sorted(pas_genes):
        s, e, strand, name = genes[lt]
        lo, hi = (s - 200, e) if strand == "+" else (s, e + 200)     # include 200 bp upstream
        o.write(f"{CHROM}\t{lo - 1}\t{hi}\t{name}\t{lt}\t{strand}\n")

with open(os.path.join(who, "analysis_key_positions.bed"), "w") as o:
    o.writelines(open(os.path.join(who, "graded_positions.bed")))
    for p, name in sorted(key):
        o.write(f"{CHROM}\t{p - 1}\t{p}\t{name}_tbdb_{p}\t{PAS}\n")
print(f"PAS: {len(rows)} literature mutations in {len(pas_genes)} genes ({', '.join(genes[g][3] for g in sorted(pas_genes))})")
