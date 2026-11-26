// Step 03: BWA-MEM to H37Rv, fixmate, sort, mark duplicates
process ALIGN {
    tag "${meta.id}"
    publishDir path: { "${params.outdir}/${meta.run_id}/${meta.id}/03_alignment" }, mode: 'copy',
               pattern: '{data,logs}/*'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("data/${meta.id}.bam"), path("data/${meta.id}.bam.bai"), emit: bam
    tuple val(meta), path('data/flagstat.json'), path('data/markdup.json'), path('logs/align_command.txt'),
          path('logs/versions.txt'), emit: stats
    path 'data/samtools_stats.txt'
    path 'logs/*'

    script:
    def ref = "${params.resources}/reference/H37Rv_NC_000962.3/H37Rv.fasta"
    """
    mkdir -p data logs
    cmd="bwa mem -t ${task.cpus} -R '@RG\\\\tID:${meta.id}\\\\tSM:${meta.id}\\\\tPL:ILLUMINA' ${ref} ${r1} ${r2}
 | samtools fixmate -m - - | samtools sort -@ 4 -o sorted.bam
 && samtools markdup -s --json -f data/markdup.json sorted.bam data/${meta.id}.bam"
    echo "\$cmd" > logs/align_command.txt
    bwa mem -t ${task.cpus} -R '@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA' ${ref} ${r1} ${r2} 2> logs/bwa.log \\
        | samtools fixmate -m - - | samtools sort -@ 4 -o sorted.bam - 2> logs/sort.log
    samtools markdup -s --json -f data/markdup.json sorted.bam data/${meta.id}.bam 2> logs/markdup.log
    samtools index data/${meta.id}.bam
    samtools flagstat -O json data/${meta.id}.bam > data/flagstat.json
    samtools stats data/${meta.id}.bam > data/samtools_stats.txt
    rm sorted.bam
    { bwa 2>&1 | grep Version; samtools --version | head -1; } > logs/versions.txt
    """
}

process ALIGN_SUMMARY {
    tag "${meta.id}"
    publishDir path: { "${params.outdir}/${meta.run_id}/${meta.id}/03_alignment" }, mode: 'copy'

    input:
    tuple val(meta), path(flagstat), path(markdup), path(cmd), path(versions)
    path thresholds

    output:
    tuple val(meta), path('03_alignment.summary.json'), emit: summary

    script:
    """
    summarise_alignment.py --flagstat ${flagstat} --markdup ${markdup} --command_file ${cmd} \\
        --versions_file ${versions} --container '${params.containers.bwa_samtools}' --sample_id ${meta.id} \\
        --run_id ${meta.run_id} --thresholds ${thresholds} --version ${workflow.manifest.version} \\
        --out 03_alignment.summary.json
    """
}
