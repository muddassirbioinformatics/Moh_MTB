// Step 12: clinical report + merged all-steps PDF, built only from the step summaries and step 09 tables
process CLINICAL_REPORT {
    tag "${meta.id}"
    publishDir path: { "${params.outdir}/${meta.run_id}/${meta.id}/12_report" }, mode: 'copy'

    input:
    tuple val(meta), path(summaries, stageAs: 'summaries/*'), path(drug_results), path(evidence),
          path(variants), path(step_pdfs, stageAs: 'steps/*')
    path drugs

    output:
    tuple val(meta), path("${meta.id}.clinical_report.pdf"), path("${meta.id}.all_steps.pdf"), emit: pdf
    tuple val(meta), path('12_report.summary.json'), emit: summary

    script:
    def mj = groovy.json.JsonOutput.toJson(meta)
    """
    cat > meta.json <<'JSON'
${mj}
JSON
    clinical_report.py --summaries summaries --drug_results ${drug_results} --evidence ${evidence} \\
        --variants ${variants} --regions ${params.resources}/who_catalogue/2023_v2/derived/analysis_regions.bed \\
        --drugs ${drugs} --meta meta.json --step_pdfs steps/*.pdf --lab_name '${params.lab_name}' \\
        --operator '${params.operator}' --banner '${params.report_banner}' --version ${workflow.manifest.version} \\
        --out_dir .
    """
}
