#!/usr/bin/env python3
"""Step 00: validate one sample's inputs and the frozen resources; write the step summary JSON.

Usage: input_check.py --sample_id ID --run_id RUN --r1 R1.fq.gz --r2 R2.fq.gz --resources DIR
                      --thresholds thresholds.yaml --version V --out 00_input.summary.json
                      [--full_resource_check]
Exit code is 0 even on FAIL; the pipeline reads "status" and stops the sample itself.
"""
import argparse, datetime, gzip, hashlib, json, os, re, subprocess
import yaml


def scan_fastq(path):
    """Decompress fully (tests gzip integrity), count reads, md5 the compressed file, keep header."""
    md5 = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            md5.update(chunk)
    n, first, maxlen, err = 0, None, 0, None
    try:
        with gzip.open(path, "rt") as fh:
            for i, line in enumerate(fh):
                if i % 4 == 0:
                    first = first or line.strip()
                elif i % 4 == 1 and n < 10000:
                    maxlen = max(maxlen, len(line.strip()))
                elif i % 4 == 3:
                    n += 1
    except (OSError, EOFError) as e:
        err = str(e)
    return {"reads": n, "md5": md5.hexdigest(), "header": first or "", "max_len": maxlen, "error": err}


def parse_header(h):
    # CASAVA 1.8+: @instrument:run:flowcell:lane:tile:x:y [read:filtered:control:index]
    f = h[1:].split()[0].split(":") if h else []
    if len(f) >= 7:
        return {"instrument": f[0], "run_number": f[1], "flowcell": f[2], "lane": f[3]}
    return None


def check_resources(res_dir, full):
    if full:
        r = subprocess.run(["sha256sum", "--quiet", "-c", "MANIFEST.sha256"], cwd=res_dir,
                           capture_output=True, text=True)
        n = sum(1 for _ in open(os.path.join(res_dir, "MANIFEST.sha256")))
        bad = [l.split(":")[0] for l in r.stdout.splitlines() if "FAILED" in l]
        return n, bad, "SHA-256 (full)"
    bad, n = [], 0
    for line in open(os.path.join(res_dir, "MANIFEST.stat")):
        p, size, mtime = line.rstrip("\n").split("\t")
        n += 1
        fp = os.path.join(res_dir, p)
        try:
            st = os.stat(fp)
            if st.st_size != int(size) or abs(st.st_mtime - float(mtime)) > 2:
                bad.append(p)
        except FileNotFoundError:
            bad.append(p + " (missing)")
    return n, bad, "size + modification time (quick)"


def main():
    ap = argparse.ArgumentParser()
    for k in ("sample_id", "run_id", "r1", "r2", "resources", "thresholds", "version", "out"):
        ap.add_argument("--" + k, required=True)
    ap.add_argument("--full_resource_check", action="store_true")
    a = ap.parse_args()
    th = yaml.safe_load(open(a.thresholds))
    t = th["00_input"]
    prov = t.get("provisional", False)
    metrics, msgs = [], []

    def add(name, plain, value, threshold, ok, unit="", warn=False):
        metrics.append({"name": name, "plain": plain, "value": value, "unit": unit, "threshold": threshold,
                        "status": "PASS" if ok else ("WARN" if warn else "FAIL"), "provisional": prov})

    add("Sample ID format", "Sample name is valid", a.sample_id, f"matches {t['sample_id_regex']}",
        bool(re.match(t["sample_id_regex"], a.sample_id)))

    r1, r2 = scan_fastq(a.r1), scan_fastq(a.r2)
    for tag, r in (("R1", r1), ("R2", r2)):
        add(f"{tag} gzip integrity", f"Sequencing file {tag[-1]} is complete and undamaged",
            "OK" if not r["error"] else "damaged", "decompresses fully", r["error"] is None)
        if r["error"]:
            msgs.append(f"{tag}: {r['error']}")
    add("R1/R2 read count equal", "Both files contain the same number of reads",
        f"{r1['reads']:,} / {r2['reads']:,}", "R1 = R2", r1["reads"] == r2["reads"])
    add("Read pairs", "Number of read pairs", r1["reads"], f">= {t['min_reads']:,}", r1["reads"] >= t["min_reads"])

    hdr = parse_header(r1["header"])
    if hdr:
        msgs.append("Run information from read header: instrument {instrument}, run {run_number}, "
                    "flow cell {flowcell}, lane {lane}.".format(**hdr))
    else:
        msgs.append(f"Read header not in Illumina format (first header: {r1['header'][:60]}); instrument, "
                    "run and flow cell not available from the file.")
    msgs.append(f"Longest read seen in first 10,000 reads: {r1['max_len']} bp.")

    n_res, bad, mode = check_resources(a.resources, a.full_resource_check)
    add("Reference resources unchanged", "Reference databases are the approved versions",
        f"{n_res - len(bad)} of {n_res} files match", f"all match MANIFEST ({mode})", not bad)
    for b in bad[:20]:
        msgs.append(f"Resource changed or missing: {b}")

    st = [m["status"] for m in metrics]
    status = "FAIL" if "FAIL" in st else "WARN" if "WARN" in st else "PASS"
    py = "Python " + ".".join(map(str, __import__("sys").version_info[:3]))
    summ = {
        "schema": "mtb-step-summary/1", "sample_id": a.sample_id, "run_id": a.run_id, "step": "00_input",
        "status": status, "pipeline_version": a.version,
        "generated": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
        "thresholds_version": th.get("version"),
        "values": {"read_pairs": f"{r1['reads']:,}", "resources_checked": n_res},
        "metrics": metrics, "messages": msgs,
        "tools": [{"name": "input_check.py", "version": a.version, "container": "mtb-report:1.1.0",
                   "command": f"input_check.py --sample_id {a.sample_id} --run_id {a.run_id} --r1 {os.path.basename(a.r1)} "
                              f"--r2 {os.path.basename(a.r2)} --resources resources/"
                              + (" --full_resource_check" if a.full_resource_check else "")},
                  {"name": "runtime", "version": py, "container": "mtb-report:1.1.0"}],
        "files": {"inputs": [{"path": os.path.basename(a.r1), "checksum": "md5:" + r1["md5"]},
                             {"path": os.path.basename(a.r2), "checksum": "md5:" + r2["md5"]}],
                  "outputs": []},
        "data": {"r1": r1, "r2": r2, "header": hdr},
    }
    json.dump(summ, open(a.out, "w"), indent=2)


if __name__ == "__main__":
    main()
