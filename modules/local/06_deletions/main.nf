// Step 06: confirm gene deletions. Runs a de novo assembly (shovill/SPAdes) only when step 04 flags
// a resistance gene with coverage drop-out, or with --always_assemble; contigs are aligned back to
// H37Rv (minimap2) to find breakpoints.
def pub06(meta) { "${params.outdir}/${meta.run_id}/${meta.id}/06_deletions" }

process SHOVILL {
    tag "${meta.id}"
    publishDir path: { pub06(meta) }, mode: 'copy', pattern: '{data,logs}/*'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path('data/contigs.fa'), emit: contigs
    path 'logs/*'

    script:
    """
    mkdir -p data logs
    shovill --R1 ${r1} --R2 ${r2} --outdir asm --gsize 4.4M --depth 100 --cpus ${task.cpus} \\
        --ram ${task.memory.toGiga() - 2} > logs/shovill.log 2>&1
    cp asm/contigs.fa data/contigs.fa
    rm -rf asm
    """
}

process CONTIG_ALIGN {
    tag "${meta.id}"
    publishDir path: { pub06(meta) }, mode: 'copy', pattern: 'data/*'

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path('data/contigs_vs_H37Rv.paf'), emit: paf

    script:
    """
    mkdir -p data
    minimap2 -cx asm5 -t ${task.cpus} ${params.resources}/reference/H37Rv_NC_000962.3/H37Rv.fasta \\
        ${contigs} > data/contigs_vs_H37Rv.paf
    """
}

process DELETION_SUMMARY {
    tag "${meta.id}"
    publishDir path: { pub06(meta) }, mode: 'copy'

    input:
    tuple val(meta), path(cov_json), path(perbase), path(asm_files)
    path thresholds

    output:
    tuple val(meta), path('06_deletions.summary.json'), emit: summary
    tuple val(meta), path('data/deletions.tsv'), emit: deletions

    script:
    def files = asm_files instanceof List ? asm_files : [asm_files]
    def paf = files.find { it.name.endsWith('.paf') }
    def ctg = files.find { it.name.endsWith('.fa') }
    """
    mkdir -p data
    find_deletions.py --coverage_summary ${cov_json} --perbase ${perbase} \\
        --regions ${params.resources}/who_catalogue/2023_v2/derived/analysis_regions.bed \\
        ${paf ? "--paf ${paf} --contigs ${ctg}" : ''} --sample_id ${meta.id} --run_id ${meta.run_id} \\
        --thresholds ${thresholds} --version ${workflow.manifest.version} \\
        --containers '${params.containers.shovill}|${params.containers.minimap2}' --out 06_deletions.summary.json
    """
}
