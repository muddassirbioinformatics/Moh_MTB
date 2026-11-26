#!/bin/bash
# Export all variants of finished samples into results/<RUN>/<SAMPLE>/variant_export/ (see tools/export_variants.py).
# Standalone: reads pipeline outputs, writes only the variant_export/ folder. Does not change the pipeline.
# Usage: tools/export_variants.sh results/<RUN>/<SAMPLE> [more sample dirs ...]
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
M=${MTBDATA_ROOT:-$HOME/mtbdata}   # symlink to the disk holding MTB_DATA (no spaces in this path)
MR=$(readlink -f "$M")
RES=$M/MTB_DATA/resources
REF=$RES/reference/H37Rv_NC_000962.3
BCF=quay.io/biocontainers/bcftools:1.24--h118bc1c_4
for d in "$@"; do
  S=$(readlink -f "$d")
  S=${S/#"$MR"/$M}
  name=$(basename "$S")
  mkdir -p "$S/variant_export"
  docker run --rm -u "$(id -u):$(id -g)" -v $M:$M $BCF bash -c \
    "bcftools csq -C 11 -p a -f $REF/H37Rv.fasta -g $REF/derived/H37Rv.csq.gff3 $S/05_variants/data/bcftools.norm.vcf.gz 2>/dev/null \
     | bcftools query -f '%POS\t%REF\t%ALT\t%QUAL\t%INFO/DP\t%INFO/DP4\t%INFO/BCSQ\n' > $S/variant_export/.genome_csq.tsv"
  python3 "$HERE/export_variants.py" "$S" "$RES" "$S/variant_export/.genome_csq.tsv"
  rm -f "$S/variant_export/.genome_csq.tsv"
done
