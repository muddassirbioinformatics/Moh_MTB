// Step 05: variant calling. bcftools genome-wide (fixed calls); LoFreq + freebayes over the resistance
// regions for minority variants; all normalised (left-aligned, split, atomised) before merging.
def pub05(meta) { "${params.outdir}/${meta.run_id}/${meta.id}/05_variants" }
def REF() { "${params.resources}/reference/H37Rv_NC_000962.3/H37Rv.fasta" }
def REGIONS() { "${params.resources}/who_catalogue/2023_v2/derived/analysis_regions.bed" }

process BCFTOOLS_CALL {
    tag "${meta.id}"
    publishDir path: { pub05(meta) }, mode: 'copy', pattern: '{data,logs}/*'

    input:
    tuple val(meta), path(bam), path(bai)
    val t

    output:
    tuple val(meta), path('raw_bcftools.vcf.gz'), emit: vcf
    path 'logs/*'

    script:
    """
    mkdir -p logs
    bcftools mpileup -f ${REF()} -q ${t.min_mapping_quality} -Q ${t.min_base_quality} -d 10000 \\
        -a AD,ADF,ADR,DP -Ou ${bam} 2> logs/mpileup.log \\
      | bcftools call -m --ploidy 1 -v -Oz -o raw_bcftools.vcf.gz 2> logs/call.log
    bcftools --version | head -1 > logs/bcftools_version.txt
    """
}

process LOFREQ {
    tag "${meta.id}"
    publishDir path: { pub05(meta) }, mode: 'copy', pattern: 'logs/*'

    input:
    tuple val(meta), path(bam), path(bai)
    val t

    output:
    tuple val(meta), path('raw_lofreq.vcf'), emit: vcf
    path 'logs/*'

    script:
    """
    mkdir -p logs
    lofreq indelqual --dindel -f ${REF()} -o iq.bam ${bam} 2> logs/indelqual.log
    lofreq index iq.bam
    lofreq call-parallel --pp-threads ${task.cpus} -f ${REF()} -l ${REGIONS()} --call-indels \\
        -q ${t.min_base_quality} -Q ${t.min_base_quality} -m ${t.min_mapping_quality} \\
        -o raw_lofreq.vcf iq.bam 2> logs/lofreq.log
    rm iq.bam*
    lofreq version | head -1 > logs/lofreq_version.txt
    """
}

process FREEBAYES {
    tag "${meta.id}"
    publishDir path: { pub05(meta) }, mode: 'copy', pattern: 'logs/*'

    input:
    tuple val(meta), path(bam), path(bai)
    val t

    output:
    tuple val(meta), path('raw_freebayes.vcf'), emit: vcf
    path 'logs/*'

    script:
    """
    mkdir -p logs
    freebayes -f ${REF()} -t ${REGIONS()} -p 1 --pooled-continuous -F 0.05 -C 3 \\
        -q ${t.min_base_quality} -m ${t.min_mapping_quality} --min-coverage ${t.min_depth} \\
        ${bam} > raw_freebayes.vcf 2> logs/freebayes.log
    freebayes --version > logs/freebayes_version.txt 2>&1
    """
}

process NORMALISE {
    tag "${meta.id}"
    publishDir path: { pub05(meta) }, mode: 'copy', pattern: 'data/*'

    input:
    tuple val(meta), path(bcf), path(lofreq), path(fb)

    output:
    tuple val(meta), path('data/bcftools.norm.vcf.gz'), path('data/lofreq.norm.vcf.gz'),
          path('data/freebayes.norm.vcf.gz'), emit: vcfs
    tuple val(meta), path('data/bcftools.norm.vcf.gz'), path('data/bcftools.norm.vcf.gz.csi'), emit: genome

    script:
    """
    mkdir -p data
    for c in bcftools:${bcf} lofreq:${lofreq} freebayes:${fb}; do
        bcftools norm -f ${REF()} -m -any -a --atom-overlaps . -Ou \${c#*:} \\
          | bcftools sort -Oz -o data/\${c%%:*}.norm.vcf.gz
        bcftools index data/\${c%%:*}.norm.vcf.gz
    done
    """
}

process VARIANT_MERGE {
    tag "${meta.id}"
    publishDir path: { pub05(meta) }, mode: 'copy'

    input:
    tuple val(meta), path(bcf), path(lofreq), path(fb)
    path thresholds

    output:
    tuple val(meta), path('05_variants.summary.json'), path('af_spectrum.png'), emit: summary
    tuple val(meta), path('data/variants.tsv'), path('data/fixed.vcf'), path('data/minority.vcf'), emit: variants

    script:
    """
    mkdir -p data
    merge_variants.py --bcftools ${bcf} --lofreq ${lofreq} --freebayes ${fb} \\
        --regions ${REGIONS()} --sample_id ${meta.id} --run_id ${meta.run_id} \\
        --thresholds ${thresholds} --version ${workflow.manifest.version} \\
        --containers '${params.containers.bcftools}|${params.containers.lofreq}|${params.containers.freebayes}' \\
        --out 05_variants.summary.json
    """
}
