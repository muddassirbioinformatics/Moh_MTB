#!/usr/bin/env bash
# Prove that a restored copy of the pipeline reproduces the original results.
# Downloads the public test isolate ERR11773715 (ENA, checks MD5), runs the pipeline into a scratch folder
# and compares with tests/expected/TEST_RUN2/TEST_BPAL_TRIAL (step summaries + result tables of the original run).
# Needs: Docker images loaded, resources in place (all databases incl. Kraken2/Sylph), source conf/env.sh done.
# Takes ~5-10 min of compute on 16 cores. Usage: setup/verify_reproduction.sh [scratch_dir]
set -euo pipefail
cd "$(dirname "$0")/.."
SCR="${1:-$PWD/verify_scratch}"; mkdir -p "$SCR/in" "$SCR/out" "$SCR/work"
ACC=ERR11773715
EXPECT_MD5_1=d85f021a50ea95ea5d7e87f269281678; EXPECT_MD5_2=f0f0b413aeff1d1fec07f3f08ad3921d
source conf/env.sh; bin/check_volume.sh
for i in 1 2; do
  f="$SCR/in/${ACC}_$i.fastq.gz"
  [[ -s "$f" ]] || wget -q -O "$f" "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR117/015/$ACC/${ACC}_$i.fastq.gz"
  exp=EXPECT_MD5_$i; [[ "$(md5sum "$f" | cut -d' ' -f1)" == "${!exp}" ]] || { echo "MD5 mismatch for $f (ENA file changed?)"; exit 1; }
done
printf 'sample_id,fastq_1,fastq_2,run_id,sample_type,sample_source\nTEST_BPAL_TRIAL,%s,%s,TEST_RUN2,gDNA,Public data ENA %s\n' \
  "$SCR/in/${ACC}_1.fastq.gz" "$SCR/in/${ACC}_2.fastq.gz" "$ACC" > "$SCR/samplesheet.csv"
nextflow run main.nf --samplesheet "$SCR/samplesheet.csv" --outdir "$SCR/out" -w "$SCR/work" -ansi-log false
python3 -I bin/compare_runs.py tests/expected/TEST_RUN2/TEST_BPAL_TRIAL "$SCR/out/TEST_RUN2/TEST_BPAL_TRIAL" \
  && echo "REPRODUCTION OK: results match the original run." || { echo "REPRODUCTION DIFFERS: see DIFF lines above."; exit 1; }
