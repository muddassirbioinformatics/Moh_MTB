// R2 (on request, --relatedness): masked consensus genomes -> pairwise SNP distances (snp-dists) and,
// with >= 4 samples, a maximum-likelihood tree (IQ-TREE)
process CONSENSUS {
    tag "${meta.id}"

    input:
    tuple val(meta), path(vcf), path(csi), path(perbase)

    output:
    path "${meta.id}.consensus.fa", emit: fasta

    script:
    def r = "${params.resources}/reference/H37Rv_NC_000962.3"
    """
    zcat ${perbase} | awk -v OFS='\\t' '\$4 < ${params.relatedness_mask_depth} {print \$1, \$2, \$3}' > low.bed
    cat ${r}/derived/low_confidence_mask.bed low.bed | sort -k1,1 -k2,2n > nmask.bed
    bcftools view -v snps -i 'GT="alt"' ${vcf} -Oz -o snps.vcf.gz && bcftools index snps.vcf.gz
    bcftools consensus -f ${r}/H37Rv.fasta -m nmask.bed snps.vcf.gz 2> consensus.log | sed "1s/.*/>${meta.id}/" > ${meta.id}.consensus.fa
    """
}

process SNP_DISTS {
    publishDir path: { "${params.outdir}/${run_id}/R2_relatedness/data" }, mode: 'copy'

    input:
    val run_id
    path fastas

    output:
    tuple path('alignment.fa'), path('snp_matrix.tsv'), emit: out

    script:
    """
    cat ${fastas} > alignment.fa
    snp-dists -b alignment.fa > snp_matrix.tsv
    """
}

process IQTREE {
    publishDir path: { "${params.outdir}/${run_id}/R2_relatedness/data" }, mode: 'copy'

    input:
    val run_id
    tuple path(aln), path(matrix)

    output:
    path 'tree.nwk'

    script:
    """
    iqtree3 -s ${aln} -m GTR+G -fast -T ${task.cpus} --prefix tree > iqtree.log 2>&1
    cp tree.treefile tree.nwk
    """
}

process RELATEDNESS_REPORT {
    publishDir path: { "${params.outdir}/${run_id}/R2_relatedness" }, mode: 'copy'

    input:
    val run_id
    tuple path(aln), path(matrix)
    path tree
    path thresholds
    path step_texts

    output:
    path 'relatedness.pdf'
    path 'R2_relatedness.summary.json'

    script:
    def t = tree.name == 'NO_TREE' ? '' : "--tree ${tree}"
    """
    relatedness.py --matrix ${matrix} ${t} --run_id ${run_id} --thresholds ${thresholds} \\
        --text ${step_texts}/R2_relatedness.yaml --version ${workflow.manifest.version} \\
        --containers '${params.containers.snp_dists}|${params.containers.iqtree}' --out_pdf relatedness.pdf
    """
}
