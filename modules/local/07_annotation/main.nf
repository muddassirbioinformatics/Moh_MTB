// Step 07: consequence annotation. bcftools csq (codon-aware for fixed calls, per-record for minority
// calls) drives loss-of-function detection; SnpEff (local H37Rv DB) supplies display names only.
def pub07(meta) { "${params.outdir}/${meta.run_id}/${meta.id}/07_annotation" }

process CSQ {
    tag "${meta.id}"
    publishDir path: { pub07(meta) }, mode: 'copy', pattern: 'data/*'

    input:
    tuple val(meta), path(tsv), path(fixed), path(minority)

    output:
    tuple val(meta), path('data/fixed.csq.vcf'), path('data/minority.csq.vcf'), emit: vcf

    script:
    def r = "${params.resources}/reference/H37Rv_NC_000962.3"
    """
    mkdir -p data
    bcftools csq -C 11 -p a -f ${r}/H37Rv.fasta -g ${r}/derived/H37Rv.csq.gff3 ${fixed} -o data/fixed.csq.vcf 2> csq_fixed.log
    bcftools csq -C 11 -l -p a -f ${r}/H37Rv.fasta -g ${r}/derived/H37Rv.csq.gff3 ${minority} -o data/minority.csq.vcf 2> csq_minority.log
    """
}

process SNPEFF {
    tag "${meta.id}"
    publishDir path: { pub07(meta) }, mode: 'copy', pattern: 'data/*'

    input:
    tuple val(meta), path(tsv), path(fixed), path(minority)

    output:
    tuple val(meta), path('data/fixed.snpeff.vcf'), path('data/minority.snpeff.vcf'), emit: vcf

    script:
    def d = "${params.resources}/snpeff/H37Rv_NC_000962.3_snpeff5.4"
    """
    mkdir -p data
    for v in fixed minority; do
        snpEff -Xmx2g -c ${d}/snpEff.config -dataDir ${d}/data -noStats -noLog -ud 200 H37Rv \${v}.vcf > data/\${v}.snpeff.vcf
    done
    """
}

process ANNOTATION_SUMMARY {
    tag "${meta.id}"
    publishDir path: { pub07(meta) }, mode: 'copy'

    input:
    tuple val(meta), path(tsv), path(fixed), path(minority), path(fcsq), path(mcsq), path(fsnp), path(msnp)
    path thresholds

    output:
    tuple val(meta), path('07_annotation.summary.json'), emit: summary
    tuple val(meta), path('data/variants.annotated.tsv'), emit: annotated

    script:
    """
    mkdir -p data
    annotate_variants.py --variants ${tsv} --csq ${fcsq} ${mcsq} --snpeff ${fsnp} ${msnp} \\
        --regions ${params.resources}/who_catalogue/2023_v2/derived/analysis_regions.bed \\
        --sample_id ${meta.id} --run_id ${meta.run_id} --thresholds ${thresholds} \\
        --version ${workflow.manifest.version} \\
        --containers '${params.containers.bcftools}|${params.containers.snpeff}' --out 07_annotation.summary.json
    """
}
