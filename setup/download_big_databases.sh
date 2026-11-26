#!/usr/bin/env bash
# Download the two large databases that were left out of the data zip, then verify against MANIFEST.sha256.
# ~8 GB (Kraken2) + ~18 GB (Sylph). Needs ~35 GB free during the Kraken2 unpack.
set -euo pipefail
R="${1:-$HOME/mtbdata/MTB_DATA/resources}"
mkdir -p "$R/kraken2/standard8_20241228" "$R/sylph/gtdb_r226_c200_dbv1"
( cd "$R/kraken2/standard8_20241228" \
  && wget -c https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20241228.tar.gz \
  && tar -xzf k2_standard_08gb_20241228.tar.gz && rm k2_standard_08gb_20241228.tar.gz ) &
( cd "$R/sylph/gtdb_r226_c200_dbv1" \
  && wget -c http://faust.compbio.cs.cmu.edu/sylph-stuff/gtdb-r226-c200-dbv1.syldb ) &
wait
echo "Verifying every database file against resources/MANIFEST.sha256 (takes a few minutes) ..."
MTBDATA_ROOT="$(dirname "$(dirname "$R")")" "$(dirname "$0")/../bin/make_manifest.sh" --check
