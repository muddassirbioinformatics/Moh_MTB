// Step 08: lineage, sub-lineage and in-silico spoligotype (TB-Profiler). TB-Profiler's own resistance
// calls are kept only as a cross-check for step 09.
def pub08(meta) { "${params.outdir}/${meta.run_id}/${meta.id}/08_lineage" }

process TBPROFILER {
    tag "${meta.id}"
    publishDir path: { pub08(meta) }, mode: 'copy', pattern: '{data,logs}/*'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path('data/tbprofiler.results.json'), emit: json
    path 'data/tbprofiler.results.txt'
    path 'logs/*'

    script:
    """
    mkdir -p data logs tmp
    tb-profiler profile -1 ${r1} -2 ${r2} --no_trim --db tbdb --db_dir ${params.resources}/tbprofiler/6.7.0 \\
        --prefix ${meta.id} --dir out --spoligotype --threads ${task.cpus} --txt --temp tmp > logs/tbprofiler.log 2>&1
    cp out/results/${meta.id}.results.json data/tbprofiler.results.json
    cp out/results/${meta.id}.results.txt data/tbprofiler.results.txt
    rm -rf out tmp
    """
}

process LINEAGE_SUMMARY {
    tag "${meta.id}"
    publishDir path: { pub08(meta) }, mode: 'copy'

    input:
    tuple val(meta), path(json)
    path thresholds

    output:
    tuple val(meta), path('08_lineage.summary.json'), emit: summary

    script:
    """
    summarise_lineage.py --tbprofiler ${json} --sample_id ${meta.id} --run_id ${meta.run_id} \\
        --thresholds ${thresholds} --version ${workflow.manifest.version} \\
        --container '${params.containers.tbprofiler}' --out 08_lineage.summary.json
    """
}
