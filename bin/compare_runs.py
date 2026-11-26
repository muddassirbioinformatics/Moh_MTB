#!/usr/bin/env python3
"""Compare two pipeline result trees (or one tree against tests/expected) for the same sample.

Usage: compare_runs.py EXPECTED_SAMPLE_DIR NEW_SAMPLE_DIR
  each argument is a folder like results/<run_id>/<sample_id> (or tests/expected/<run_id>/<sample_id>).

What is compared (only files present in EXPECTED):
  * every <step>.summary.json: all fields EXCEPT "generated" (timestamp), and with output-folder paths removed
  * data tables: 09_resistance/data/drug_results.tsv + evidence.tsv, 05_variants/data/variants.tsv,
    07_annotation/data/variants.annotated.tsv (byte-for-byte after line sort)
PDFs and logs are not compared (they carry timestamps and paths). Exit code 0 = identical, 1 = differences.
"""
import json, os, re, sys

TABLES = ["09_resistance/data/drug_results.tsv", "09_resistance/data/evidence.tsv",
          "05_variants/data/variants.tsv", "07_annotation/data/variants.annotated.tsv"]

def norm(o, roots):
    if isinstance(o, dict):
        return {k: norm(v, roots) for k, v in o.items() if k != "generated"}
    if isinstance(o, list):
        return [norm(v, roots) for v in o]
    if isinstance(o, str):
        for r in roots:
            o = o.replace(r, "<OUT>")
        o = re.sub(r"/[^ \"']*/work/[0-9a-f]{2}/[0-9a-f]{30}", "<WORK>", o)
        return re.sub(r"[^ \"']*/(MTB_DATA|MTB_PIPELINE)/", "<ROOT>/", o)
    return o

def diff(a, b, path=""):
    out = []
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            if k not in a: out.append(f"{path}/{k}: only in NEW")
            elif k not in b: out.append(f"{path}/{k}: only in EXPECTED")
            else: out += diff(a[k], b[k], f"{path}/{k}")
    elif isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
        for i, (x, y) in enumerate(zip(a, b)): out += diff(x, y, f"{path}[{i}]")
    elif a != b:
        out.append(f"{path}: expected {str(a)[:80]!r} got {str(b)[:80]!r}")
    return out

def main(exp, new):
    exp, new = exp.rstrip("/"), new.rstrip("/")
    roots = [os.path.dirname(os.path.dirname(exp)), os.path.dirname(os.path.dirname(new))]
    n_ok, problems = 0, []
    for d, _, fs in os.walk(exp):
        for f in sorted(fs):
            rel = os.path.relpath(os.path.join(d, f), exp)
            if f.endswith(".summary.json"):
                p = os.path.join(new, rel)
                if not os.path.exists(p): problems.append(f"{rel}: missing in NEW"); continue
                a = norm(json.load(open(os.path.join(d, f))), roots)
                b = norm(json.load(open(p)), roots)
                dd = diff(a, b, rel)
                problems += dd; n_ok += not dd
            elif rel in TABLES:
                p = os.path.join(new, rel)
                if not os.path.exists(p): problems.append(f"{rel}: missing in NEW"); continue
                a = sorted(open(os.path.join(d, f)).read().splitlines())
                b = sorted(open(p).read().splitlines())
                if a != b: problems.append(f"{rel}: {len(set(a) ^ set(b))} differing lines")
                else: n_ok += 1
    print(f"identical files: {n_ok}; differences: {len(problems)}")
    for p in problems[:60]: print("  DIFF", p)
    return 1 if problems else 0

if __name__ == "__main__":
    if len(sys.argv) != 3: sys.exit(__doc__)
    sys.exit(main(sys.argv[1], sys.argv[2]))
