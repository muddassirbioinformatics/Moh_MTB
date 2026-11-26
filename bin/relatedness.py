#!/usr/bin/env python3
"""R2: pairwise SNP distances -> transmission categories (Walker et al. 2013) and a run-level PDF."""
import argparse, csv, datetime, itertools
import yaml
from mtbstep import Step
from step_report import build


def clusters(names, d, cut):
    """Single-linkage groups of samples within `cut` SNPs."""
    parent = {n: n for n in names}

    def find(x):
        while parent[x] != x:
            x = parent[x]
        return x
    for a, b in itertools.combinations(names, 2):
        if d[a][b] <= cut:
            parent[find(a)] = find(b)
    groups = {}
    for n in names:
        groups.setdefault(find(n), []).append(n)
    return [g for g in groups.values() if len(g) > 1]


def main():
    ap = argparse.ArgumentParser()
    for k in ("matrix", "run_id", "thresholds", "text", "version", "containers", "out_pdf"):
        ap.add_argument("--" + k, required=True)
    ap.add_argument("--tree")
    a = ap.parse_args()
    s = Step("R2_relatedness", f"run {a.run_id}", a.run_id, a.version, a.thresholds)
    t = s.t
    rows = list(csv.reader(open(a.matrix), delimiter="\t"))
    names = rows[0][1:]
    d = {r[0]: dict(zip(names, map(int, r[1:]))) for r in rows[1:]}
    pairs = []
    for x, y in itertools.combinations(names, 2):
        n = d[x][y]
        cat = ("likely recent transmission" if n <= t["linked_max_snps"] else
               "possible" if n < t["unlinked_min_snps"] else "unlikely")
        pairs.append([x, y, str(n), cat])
    pairs.sort(key=lambda p: int(p[2]))
    linked = clusters(names, d, t["linked_max_snps"])
    s.d["results"] = {"title": "Closest pairs", "columns": ["Sample A", "Sample B", "SNPs apart", "Recent transmission"],
                      "rows": pairs[:25], "widths_mm": [45, 45, 30, 60]}
    s.d["tables"] = [{"title": "Pairwise SNP distance matrix", "columns": [""] + names,
                      "rows": [[n] + [str(d[n][m]) for m in names] for n in names]}]
    s.metric("Samples compared", "Samples compared", len(names))
    s.metric("Linked clusters", "Groups of possibly linked patients",
             "; ".join(", ".join(g) for g in linked) or "none",
             threshold=f"<= {t['linked_max_snps']} SNPs", status="WARN" if linked else "PASS")
    s.msg(f"Categories (Walker et al. 2013): <= {t['linked_max_snps']} SNPs likely recent transmission; "
          f"{t['linked_max_snps'] + 1}-{t['unlinked_min_snps'] - 1} possible; >= {t['unlinked_min_snps']} unlikely.")
    s.msg("Consensus genomes exclude the low-confidence mask (repeats, PE/PPE) and positions below "
          f"{t['mask_depth_x']}x; only SNPs are compared. Genomic links support, but do not replace, "
          "epidemiological investigation.")
    if a.tree:
        s.msg("Maximum-likelihood tree (IQ-TREE, GTR+G): " + open(a.tree).read().strip()[:600])
    c_s, c_i = a.containers.split("|")
    s.tool("snp-dists", c_s.split(":")[-1].split("--")[0], c_s, "snp-dists -b alignment.fa")
    if a.tree:
        s.tool("IQ-TREE", c_i.split(":")[-1].split("--")[0], c_i, "iqtree3 -m GTR+G -fast")
    s.write("R2_relatedness.summary.json", n=len(names), clusters="; ".join(", ".join(g) for g in linked) or "none")
    s.d["generated"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    build(s.d, yaml.safe_load(open(a.text)), a.out_pdf, [])


if __name__ == "__main__":
    main()
