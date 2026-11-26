// Step 02: species identification (Mykrobe + Sylph must agree) and contamination screen (Kraken2)
def pub02(meta) { "${params.outdir}/${meta.run_id}/${meta.id}/02_species" }

process MYKROBE {
    tag "${meta.id}"
    publishDir path: { pub02(meta) }, mode: 'copy', pattern: '{data,logs}/*'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path('data/mykrobe.json'), emit: json
    path 'logs/*'

    script:
    """
    mkdir -p data logs tmp
    mykrobe predict -s ${meta.id} -S tb --panels_dir ${params.resources}/mykrobe/0.13.0 \\
        -1 ${r1} ${r2} --format json -o data/mykrobe.json -t ${task.cpus} --tmp tmp/ 2> logs/mykrobe.log
    mykrobe --version > logs/mykrobe_version.txt 2>&1 || true
    """
}

process SYLPH {
    tag "${meta.id}"
    publishDir path: { pub02(meta) }, mode: 'copy', pattern: '{data,logs}/*'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path('data/sylph.tsv'), emit: tsv
    path 'logs/*'

    script:
    """
    mkdir -p data logs
    sylph profile ${params.resources}/sylph/gtdb_r226_c200_dbv1/gtdb-r226-c200-dbv1.syldb \\
        -1 ${r1} -2 ${r2} -t ${task.cpus} -o data/sylph.tsv 2> logs/sylph.log
    """
}

process KRAKEN2 {
    tag "${meta.id}"
    publishDir path: { pub02(meta) }, mode: 'copy', pattern: '{data,logs}/*'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path('data/kraken2.report'), emit: report
    path 'logs/*'

    script:
    """
    mkdir -p data logs
    kraken2 --db ${params.resources}/kraken2/standard8_20241228 --paired --gzip-compressed \\
        --threads ${task.cpus} --report data/kraken2.report --output /dev/null ${r1} ${r2} 2> logs/kraken2.log
    """
}

process SPECIES_SUMMARY {
    tag "${meta.id}"
    publishDir path: { pub02(meta) }, mode: 'copy'

    input:
    tuple val(meta), path(myk), path(syl), path(k2)
    path thresholds

    output:
    tuple val(meta), path('02_species.summary.json'), emit: summary

    script:
    """
    summarise_species.py --mykrobe ${myk} --sylph ${syl} --kraken2 ${k2} --sample_id ${meta.id} \\
        --run_id ${meta.run_id} --thresholds ${thresholds} --version ${workflow.manifest.version} \\
        --containers '${params.containers.mykrobe}|${params.containers.sylph}|${params.containers.kraken2}' \\
        --out 02_species.summary.json
    """
}
