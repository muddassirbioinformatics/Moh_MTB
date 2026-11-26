// Step 01: read QC and trimming with fastp; summary JSON + quality plot
process FASTP {
    tag "${meta.id}"
    publishDir path: { "${params.outdir}/${meta.run_id}/${meta.id}/01_read_qc" }, mode: 'copy', pattern: '{data,logs}/*'

    input:
    tuple val(meta), path(r1), path(r2)
    val fp      // fastp settings from conf/thresholds.yaml

    output:
    tuple val(meta), path("${meta.id}_R1.trim.fastq.gz"), path("${meta.id}_R2.trim.fastq.gz"), emit: reads
    tuple val(meta), path('data/fastp.json'), path('logs/fastp_command.txt'), path('logs/fastp_version.txt'), emit: report
    path 'data/fastp.html'
    path 'logs/*'

    script:
    """
    mkdir -p data logs
    cmd="fastp -i ${r1} -I ${r2} -o ${meta.id}_R1.trim.fastq.gz -O ${meta.id}_R2.trim.fastq.gz \\
 --detect_adapter_for_pe --cut_right --cut_right_window_size ${fp.cut_right_window_size} \\
 --cut_right_mean_quality ${fp.cut_right_mean_quality} --length_required ${fp.length_required} \\
 --trim_poly_g -w ${task.cpus} -j data/fastp.json -h data/fastp.html"
    echo "\$cmd" > logs/fastp_command.txt
    eval "\$cmd" 2> logs/fastp.log
    fastp --version > logs/fastp_version.txt 2>&1
    cp .command.sh logs/command.sh
    """
}

process READ_QC_SUMMARY {
    tag "${meta.id}"
    publishDir path: { "${params.outdir}/${meta.run_id}/${meta.id}/01_read_qc" }, mode: 'copy'

    input:
    tuple val(meta), path(json), path(cmd), path(ver)
    path thresholds

    output:
    tuple val(meta), path('01_read_qc.summary.json'), path('quality.png'), emit: summary

    script:
    """
    read_qc_summary.py --fastp ${json} --sample_id ${meta.id} --run_id ${meta.run_id} \\
        --thresholds ${thresholds} --version ${workflow.manifest.version} \\
        --fastp_version_file ${ver} --fastp_container '${params.containers.fastp}' \\
        --command_file ${cmd} --out 01_read_qc.summary.json --plot quality.png
    """
}
