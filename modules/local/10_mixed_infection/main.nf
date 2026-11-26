// Step 10: mixed infection. Genome-wide heterozygous-site test outside the low-confidence mask
// (MixInfect-style), lineage mixture from TB-Profiler, and the OmegaTB rule for comparison.
// The mask is applied in mixed_infection.py: bcftools 1.24 mpileup segfaults with -T ^mask.bed.
process HET_SITES {
    tag "${meta.id}"
    publishDir path: { "${params.outdir}/${meta.run_id}/${meta.id}/10_mixed_infection" }, mode: 'copy',
               pattern: 'data/*'

    input:
    tuple val(meta), path(bam), path(bai)
    val t

    output:
    tuple val(meta), path('data/het_q20.vcf.gz'), path('data/het_omega.vcf.gz'), emit: vcf

    script:
    def ref = "${params.resources}/reference/H37Rv_NC_000962.3/H37Rv.fasta"
    """
    mkdir -p data
    for q in 20:data/het_q20.vcf.gz ${t.omega_min_mapq}:data/het_omega.vcf.gz; do
        bcftools mpileup -f ${ref} -q \${q%%:*} -Q 20 -d 10000 -a AD,ADF,ADR,DP -Ou ${bam} \\
          | bcftools call -m -v -V indels -Oz -o \${q#*:}
    done
    """
}

process MIXED_SUMMARY {
    tag "${meta.id}"
    publishDir path: { "${params.outdir}/${meta.run_id}/${meta.id}/10_mixed_infection" }, mode: 'copy'

    input:
    tuple val(meta), path(q20), path(omega), path(tbprofiler)
    path thresholds

    output:
    tuple val(meta), path('10_mixed_infection.summary.json'), path('minor_af.png'), emit: summary

    script:
    """
    mixed_infection.py --mask ${params.resources}/reference/H37Rv_NC_000962.3/derived/low_confidence_mask.bed \
        --vcf ${q20} --omega_vcf ${omega} --tbprofiler ${tbprofiler} --sample_id ${meta.id} \\
        --run_id ${meta.run_id} --thresholds ${thresholds} --version ${workflow.manifest.version} \\
        --container '${params.containers.bcftools}' --out 10_mixed_infection.summary.json
    """
}
