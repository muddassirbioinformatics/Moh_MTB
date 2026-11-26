// Step 00: validate sample inputs and frozen resources (summary JSON only; PDF via STEP_REPORT)
process INPUT_CHECK {
    tag "${meta.id}"
    publishDir path: { "${params.outdir}/${meta.run_id}/${meta.id}/00_input" }, mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2)
    path thresholds

    output:
    tuple val(meta), path('00_input.summary.json'), emit: summary
    path 'logs/*'

    script:
    def full = params.full_resource_check ? '--full_resource_check' : ''
    """
    mkdir -p logs
    input_check.py --sample_id ${meta.id} --run_id ${meta.run_id} --r1 ${r1} --r2 ${r2} \\
        --resources ${params.resources} --thresholds ${thresholds} \\
        --version ${workflow.manifest.version} --out 00_input.summary.json ${full} 2> logs/input_check.log
    cp .command.sh logs/command.sh
    """
}
