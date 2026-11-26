// Step 09: per-drug resistance grading against the WHO 2023 catalogue
process GRADE_RESISTANCE {
    tag "${meta.id}"
    publishDir path: { "${params.outdir}/${meta.run_id}/${meta.id}/09_resistance" }, mode: 'copy'

    input:
    tuple val(meta), path(variants), path(deletions), path(cov_json), path(perbase), path(tbprofiler)
    path drugs
    path thresholds

    output:
    tuple val(meta), path('09_resistance.summary.json'), emit: summary
    tuple val(meta), path('data/drug_results.tsv'), path('data/evidence.tsv'), emit: results

    script:
    """
    mkdir -p data
    grade_resistance.py --variants ${variants} --deletions ${deletions} --coverage_summary ${cov_json} \\
        --perbase ${perbase} --derived ${params.resources}/who_catalogue/2023_v2/derived --drugs ${drugs} \\
        --tbprofiler ${tbprofiler} --sample_id ${meta.id} --run_id ${meta.run_id} --thresholds ${thresholds} \\
        --version ${workflow.manifest.version} --out 09_resistance.summary.json
    """
}
