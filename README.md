# MTB_PIPELINE: whole-genome TB drug-resistance pipeline

A Nextflow pipeline that takes paired-end Illumina FASTQ from a cultured *Mycobacterium tuberculosis*
complex (MTBC) isolate and produces a whole-genome drug-resistance report for 17 drugs, graded against
the WHO 2023 mutation catalogue (v2), plus species check, lineage/spoligotype, mixed-infection check,
WHO category (Hr-TB / RR-TB / MDR / pre-XDR / XDR) and optional transmission relatedness.

**Status: validation build (0.1.0-dev).** Every clinical report carries the banner
"VALIDATION BUILD - NOT FOR CLINICAL USE" until the pipeline's own validation (V1-V6) is signed off.

## Run
    source conf/env.sh                       # Java 17 + Nextflow 26.10.4 (pinned)
    bin/check_volume.sh                      # checks the data disk against conf/base.config
    nextflow run main.nf --samplesheet raw/<RUN>/samplesheet.csv --lab_name "<Lab name>" [--relatedness]

Samplesheet (CSV): `sample_id,fastq_1,fastq_2,run_id` plus optional lab columns printed on the report:
`client_id, test_id, patient_name, age, sex, referring_doctor, sample_type, sample_source, receipt_date,
sample_temperature, sample_qc_date, dna_260_280, dna_conc_ng_ul, library_type, library_date, library_qc,
sequencer, run_date`. Sample IDs: letters, digits, `_` and `-` only.

Options: `--always_assemble` (step 06 assembly for every sample), `--full_resource_check` (SHA-256 of all
databases instead of the quick size/date check), `--relatedness` (R2), `--operator`, `--report_banner`.

## Steps and outputs
`results/<run_id>/<sample_id>/<NN_step>/` holds `<step>.summary.json`, `<sample>.<step>.pdf`
(plain-language page + technical pages), `data/` and `logs/`.

| Step | What | Main tools |
|---|---|---|
| 00_input | files, naming, database integrity | input_check.py |
| 01_read_qc | trimming, yield, Q30 | fastp |
| 02_species | MTBC (two methods), other organisms | Mykrobe, Sylph, Kraken2 (screen) |
| 03_alignment | map to H37Rv | BWA-MEM, samtools |
| 04_coverage | depth, breadth, gene drop-out, key positions | mosdepth |
| 05_variants | fixed + minority variants | bcftools, LoFreq, freebayes |
| 06_deletions | deletion confirmation (only on drop-out) | shovill, minimap2 |
| 07_annotation | consequences, LoF | bcftools csq, SnpEff |
| 08_lineage | lineage, spoligotype | TB-Profiler |
| 09_resistance | WHO 2023 grading, 17 drugs | grade_resistance.py |
| 10_mixed_infection | mixed strains | bcftools + mixed_infection.py |
| 11_classification | Hr/RR/MDR/pre-XDR/XDR | classify.py |
| 12_report | clinical report + merged all-steps PDF | clinical_report.py |
| R1_run_qc | run summary + MultiQC | run_summary.py, MultiQC |
| R2_relatedness | SNP distances, tree (on request) | snp-dists, IQ-TREE |

## Setup
`resources/`, `raw/`, `results/` and `validation/data` are not part of this repo (see `.gitignore`) — point
them at wherever you keep sequencing data and reference databases, for example via symlinks. See `setup/`
for installing tools (`install_tools.sh`), wiring up data paths (`setup_links.sh`), loading pinned containers
(`load_containers.sh`), downloading the large databases (`download_big_databases.sh`), and checking a fresh
checkout reproduces the bundled `tests/expected/` output (`verify_reproduction.sh`).

Pinned containers: `conf/containers.tsv`. Database versions: `conf/resources.yaml`. All cut-offs:
`conf/thresholds.yaml` (marked provisional until validation). Derived tables are rebuilt with
`bin/build_who_tables.py`, `bin/build_analysis_targets.py`, `bin/build_csq_gff.py`, then
`bin/make_manifest.sh --add <dir>`.

`tools/export_variants.sh` exports every variant of every sample (resistant, not detected, or VUS) as a
per-sample TSV, independent of the main pipeline run.

## License
MIT — see `LICENSE`.
