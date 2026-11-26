// Step 04: genome and resistance-gene coverage (mosdepth); per-base depth kept for step 09
process MOSDEPTH {
    tag "${meta.id}"
    publishDir path: { "${params.outdir}/${meta.run_id}/${meta.id}/04_coverage" }, mode: 'copy',
               pattern: '{data,logs}/*'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path('data/genome.mosdepth.summary.txt'), path('data/genome.mosdepth.global.dist.txt'),
          path('data/genome.per-base.bed.gz'), path('data/genes.regions.bed.gz'),
          path('data/genes.thresholds.bed.gz'), path('data/key.regions.bed.gz'), emit: cov
    tuple val(meta), path('data/genome.per-base.bed.gz'), emit: perbase
    path 'logs/*'

    script:
    def d = "${params.resources}/who_catalogue/2023_v2/derived"
    """
    mkdir -p data logs
    mosdepth -t ${task.cpus} data/genome ${bam} 2> logs/mosdepth_genome.log
    mosdepth -t ${task.cpus} -n --by ${d}/analysis_regions.bed -T 1,10,50 data/genes ${bam} 2> logs/mosdepth_genes.log
    mosdepth -t ${task.cpus} -n --by ${d}/analysis_key_positions.bed data/key ${bam} 2> logs/mosdepth_key.log
    mosdepth --version > logs/mosdepth_version.txt 2>&1
    """
}

process COVERAGE_SUMMARY {
    tag "${meta.id}"
    publishDir path: { "${params.outdir}/${meta.run_id}/${meta.id}/04_coverage" }, mode: 'copy'

    input:
    tuple val(meta), path(summary), path(dist), path(perbase), path(genes), path(gthr), path(key)
    path thresholds

    output:
    tuple val(meta), path('04_coverage.summary.json'), path('genome_depth.png'), path('gene_depth.png'), emit: summary

    script:
    """
    summarise_coverage.py --summary ${summary} --dist ${dist} --perbase ${perbase} --genes ${genes} \\
        --gene_thresholds ${gthr} --key ${key} --container '${params.containers.mosdepth}' \\
        --sample_id ${meta.id} --run_id ${meta.run_id} --thresholds ${thresholds} \\
        --version ${workflow.manifest.version} --out 04_coverage.summary.json
    """
}
