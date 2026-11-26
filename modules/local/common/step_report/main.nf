// Renders <sample>.<step>.pdf (plain-language page + technical pages) from a step summary JSON
process STEP_REPORT {
    tag "${meta.id}:${step}"
    publishDir path: { "${params.outdir}/${meta.run_id}/${meta.id}/${step}" }, mode: 'copy'

    input:
    tuple val(meta), val(step), path(summary), path(plots)
    path step_texts

    output:
    tuple val(meta), val(step), path("${meta.id}.${step}.pdf"), emit: pdf

    script:
    def pl = plots instanceof List ? plots : [plots]
    def p = pl ? "--plots ${pl.join(' ')}" : ''
    """
    step_report.py ${summary} ${step_texts}/${step}.yaml ${meta.id}.${step}.pdf ${p}
    """
}
