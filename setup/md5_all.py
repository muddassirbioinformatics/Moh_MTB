#!/usr/bin/env python3
"""MD5 (and SHA-256) of every file inside each zip of the export, plus of the zips themselves.
Usage: md5_all.py EXPORT_DIR      -> writes MD5SUMS.txt (zips + guide) and MD5_<zipname>_contents.txt (one line per file inside)
Verify later:  md5sum -c MD5SUMS.txt   (in EXPORT_DIR);  inside an unzipped tree:  md5sum -c <(sed 's/^\\([0-9a-f]*\\)  /\\1  /' MD5_<zip>_contents.txt)
"""
import hashlib, os, sys, zipfile

def md5_file(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""): h.update(b)
    return h.hexdigest()

d = sys.argv[1]
names = sorted(f for f in os.listdir(d) if f.endswith(".zip") or f == "REPLICATE_PIPELINE.md" or f.startswith("MTB_PIPELINE_Build_Chat_Log"))
with open(os.path.join(d, "MD5SUMS.txt"), "w") as out:
    for n in names: out.write(f"{md5_file(os.path.join(d, n))}  {n}\n")
for n in names:
    if not n.endswith(".zip"): continue
    z = zipfile.ZipFile(os.path.join(d, n)); cnt = 0
    with open(os.path.join(d, f"MD5_{n[:-4]}_contents.txt"), "w") as out:
        for i in sorted(z.infolist(), key=lambda i: i.filename):
            if i.is_dir() or (i.external_attr >> 28) == 0xA:   # skip dirs and symlinks (listed separately)
                if not i.is_dir(): out.write(f"SYMLINK  {i.filename}  -> {z.read(i).decode()}\n")
                continue
            h = hashlib.md5()
            with z.open(i) as f:
                for b in iter(lambda: f.read(1 << 22), b""): h.update(b)
            out.write(f"{h.hexdigest()}  {i.filename}\n"); cnt += 1
    print(n, cnt, "files hashed")
