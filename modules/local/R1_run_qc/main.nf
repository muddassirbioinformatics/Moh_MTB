// R1: run-level QC. MultiQC over all samples + a run summary PDF.
process MULTIQC {
    publishDir path: { "${params.outdir}/${run_id}/R1_run_qc" }, mode: 'copy'

    input:
    val run_id
    path files, stageAs: 'inputs/?/*'

    output:
    path 'multiqc.html'

    script:
    """
    multiqc -f -n multiqc.html inputs > multiqc.log 2>&1
    """
}

process RUN_SUMMARY {
    publishDir path: { "${params.outdir}/${run_id}/R1_run_qc" }, mode: 'copy'

    input:
    val run_id
    path summaries, stageAs: 'summaries/?/*'
    path step_texts

    output:
    path 'run_qc.pdf'
    path 'R1_run_qc.summary.json'

    script:
    """
    mkdir -p flat && for f in summaries/*/*.json; do cp \$f flat/\$(echo \$f | tr '/' '_'); done
    run_summary.py --summaries flat --text ${step_texts}/R1_run_qc.yaml --version ${workflow.manifest.version} \\
        --out_pdf run_qc.pdf
    """
}
