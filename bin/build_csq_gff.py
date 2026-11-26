#!/usr/bin/env python3
"""Convert the NCBI H37Rv GFF3 into the Ensembl-style layout bcftools csq expects
(gene:ID / transcript:ID / CDS+exon children). Run once per reference version.

Usage: build_csq_gff.py NCBI.gff3 OUT.gff3
"""
import sys

src, out = sys.argv[1:3]
genes, children = {}, {}
for line in open(src):
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 9:
        continue
    a = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
    if f[2] in ("gene", "pseudogene"):
        genes[a["ID"]] = (f, a)
    elif f[2] in ("CDS", "rRNA", "tRNA", "ncRNA") and a.get("Parent", "").startswith("gene-"):
        children.setdefault(a["Parent"], []).append((f, a))

n = 0
with open(out, "w") as o:
    o.write("##gff-version 3\n")
    for gid, (gf, ga) in genes.items():
        kids = children.get(gid, [])
        if not kids:
            continue
        lt = ga.get("locus_tag", gid[5:])
        name = ga.get("gene", lt)
        coding = kids[0][0][2] == "CDS"
        biotype = "protein_coding" if coding else kids[0][0][2]
        base = gf[:8]
        o.write("\t".join(base[:2] + ["gene"] + base[3:]) + f"\tID=gene:{lt};Name={name};biotype={biotype}\n")
        o.write("\t".join(base[:2] + ["transcript"] + base[3:]) +
                f"\tID=transcript:{lt};Parent=gene:{lt};Name={name};biotype={biotype}\n")
        for kf, _ in kids:
            for typ in (("CDS",) if coding else ("exon",)):
                phase = kf[7] if typ == "CDS" else "."
                o.write("\t".join(kf[:2] + [typ] + kf[3:7] + [phase]) + f"\tParent=transcript:{lt}\n")
        n += 1
print(f"wrote {n} genes")
