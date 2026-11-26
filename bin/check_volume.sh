#!/usr/bin/env bash
# Refuse to run unless ~/mtbdata is the MTB_DATA disk (UUID 8A3ABC8D3ABC77B1).
# Both external disks carry the label "New Volume", so the mount name alone is not trusted.
set -euo pipefail
EXPECTED_UUID="8A3ABC8D3ABC77B1"
ROOT="${MTBDATA_ROOT:-$HOME/mtbdata}"
src=$(findmnt -no SOURCE --target "$(readlink -f "$ROOT")") || { echo "ERROR: $ROOT is not mounted" >&2; exit 1; }
uuid=$(lsblk -no UUID "$src")
if [[ "$uuid" != "$EXPECTED_UUID" ]]; then
  echo "ERROR: $ROOT is on $src (UUID $uuid), expected UUID $EXPECTED_UUID" >&2; exit 1
fi
grep -q "UUID=$EXPECTED_UUID" "$ROOT/MTB_DATA/.volume_id" || { echo "ERROR: .volume_id marker missing" >&2; exit 1; }
echo "OK: MTB_DATA volume $src ($uuid)"
