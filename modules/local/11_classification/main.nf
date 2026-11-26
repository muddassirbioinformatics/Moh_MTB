// Step 11: WHO programmatic category (Hr / RR / MDR / pre-XDR / XDR)
process CLASSIFY {
    tag "${meta.id}"
    publishDir path: { "${params.outdir}/${meta.run_id}/${meta.id}/11_classification" }, mode: 'copy'

    input:
    tuple val(meta), path(drug_results), path(mixed)
    path thresholds

    output:
    tuple val(meta), path('11_classification.summary.json'), emit: summary

    script:
    """
    classify.py --drug_results ${drug_results} --mixed ${mixed} --sample_id ${meta.id} --run_id ${meta.run_id} \\
        --thresholds ${thresholds} --version ${workflow.manifest.version} --out 11_classification.summary.json
    """
}
