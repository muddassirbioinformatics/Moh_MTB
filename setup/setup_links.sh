#!/usr/bin/env bash
# Recreate the data symlinks of the pipeline on a new machine.
# Usage: setup/setup_links.sh /path/to/disk_or_folder      (the folder that will CONTAIN MTB_DATA/)
#   - creates <path>/MTB_DATA/{resources,raw,results,validation,archive,containers}
#   - writes the .volume_id marker (UUID taken from the mounted disk, or pass UUID=<uuid> in the environment)
#   - links ~/mtbdata -> <path>   and   ./resources ./raw ./results ./validation/data -> ~/mtbdata/MTB_DATA/...
# bin/check_volume.sh compares the disk UUID with EXPECTED_UUID; on a new disk update that value too.
set -euo pipefail
DEST="${1:?usage: setup_links.sh /path/that/will/contain/MTB_DATA}"
UUID_VAL="${UUID:-$(findmnt -no UUID --target "$DEST" 2>/dev/null || echo UNKNOWN)}"
mkdir -p "$DEST/MTB_DATA"/{resources,raw,results,validation,archive,containers}
printf 'MTB_DATA volume\nUUID=%s\ncreated=%s\n' "$UUID_VAL" "$(date -I)" > "$DEST/MTB_DATA/.volume_id"
ln -sfn "$DEST" "$HOME/mtbdata"
cd "$(dirname "$0")/.."
mkdir -p validation
for l in resources raw results; do ln -sfn "$HOME/mtbdata/MTB_DATA/$l" "$l"; done
ln -sfn "$HOME/mtbdata/MTB_DATA/validation" validation/data
echo "Linked. ~/mtbdata -> $DEST (UUID $UUID_VAL)."
echo "Now put the data zip content into $DEST/MTB_DATA/, then edit EXPECTED_UUID in bin/check_volume.sh to $UUID_VAL."
