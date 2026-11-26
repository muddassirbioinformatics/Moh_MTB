#!/usr/bin/env bash
# Write (default) or verify (--check) resources/MANIFEST.sha256 and MANIFEST.stat.
# --add DIR adds the files under DIR (relative to resources/) without re-hashing everything.
# A failed check means a database changed.
set -euo pipefail
cd "${MTBDATA_ROOT:-$HOME/mtbdata}/MTB_DATA/resources"
if [[ "${1:-}" == "--add" ]]; then
  d="./${2:?usage: --add DIR}"
  { grep -v "  ${d}/" MANIFEST.sha256 || true; find "$d" -type f -print0 | xargs -0 sha256sum; } | sort -k2 > MANIFEST.tmp
  mv MANIFEST.tmp MANIFEST.sha256
  { grep -v "^${d}/" MANIFEST.stat || true; find "$d" -type f -printf '%p\t%s\t%T@\n'; } | sort > MANIFEST.tmp
  mv MANIFEST.tmp MANIFEST.stat
  echo "added $d; manifest now $(wc -l < MANIFEST.sha256) files"
elif [[ "${1:-}" == "--check" ]]; then
  sha256sum --quiet -c MANIFEST.sha256 && echo "OK: all resource checksums match"
else
  find . -type f ! -name 'MANIFEST.*' -print0 | sort -z | xargs -0 sha256sum > MANIFEST.sha256
  # quick per-run check: path, size, mtime (step 00 compares against this)
  find . -type f ! -name 'MANIFEST.*' -printf '%p\t%s\t%T@\n' | sort > MANIFEST.stat
  echo "wrote MANIFEST.sha256 and MANIFEST.stat ($(wc -l < MANIFEST.sha256) files)"
fi
