#!/usr/bin/env nextflow
// MTB_PIPELINE: whole-genome TB drug-resistance pipeline. Steps are numbered as in
// docs/planning/Pipeline_Architecture_and_Storage_Plan.pdf; each step writes a summary JSON and a PDF.
nextflow.enable.dsl = 2

include { INPUT_CHECK }               from './modules/local/00_input/main.nf'
include { FASTP; READ_QC_SUMMARY }    from './modules/local/01_read_qc/main.nf'
include { MYKROBE; SYLPH; KRAKEN2; SPECIES_SUMMARY } from './modules/local/02_species/main.nf'
include { ALIGN; ALIGN_SUMMARY }      from './modules/local/03_alignment/main.nf'
include { MOSDEPTH; COVERAGE_SUMMARY } from './modules/local/04_coverage/main.nf'
include { BCFTOOLS_CALL; LOFREQ; FREEBAYES; NORMALISE; VARIANT_MERGE } from './modules/local/05_variants/main.nf'
include { SHOVILL; CONTIG_ALIGN; DELETION_SUMMARY } from './modules/local/06_deletions/main.nf'
include { CSQ; SNPEFF; ANNOTATION_SUMMARY } from './modules/local/07_annotation/main.nf'
include { TBPROFILER; LINEAGE_SUMMARY } from './modules/local/08_lineage/main.nf'
include { GRADE_RESISTANCE } from './modules/local/09_resistance/main.nf'
include { HET_SITES; MIXED_SUMMARY } from './modules/local/10_mixed_infection/main.nf'
include { CLASSIFY } from './modules/local/11_classification/main.nf'
include { CLINICAL_REPORT } from './modules/local/12_report/main.nf'
include { MULTIQC; RUN_SUMMARY } from './modules/local/R1_run_qc/main.nf'
include { CONSENSUS; SNP_DISTS; IQTREE; RELATEDNESS_REPORT } from './modules/local/R2_relatedness/main.nf'
include { STEP_REPORT }               from './modules/local/common/step_report/main.nf'

def status(json) { new groovy.json.JsonSlurper().parse(json.toFile()).status }
def dropouts(json) { new groovy.json.JsonSlurper().parse(json.toFile()).values.dropout_genes ?: [] }

workflow {
    if (!params.samplesheet) error "Provide --samplesheet (CSV: sample_id,fastq_1,fastq_2,run_id)"
    def TH = new org.yaml.snakeyaml.Yaml().load(file(params.thresholds).text)
    thresholds = file(params.thresholds)
    step_texts = file(params.step_texts)

    reads = channel.fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { r -> tuple([id: r.sample_id, run_id: r.run_id,
                           lims: r.findAll { k, v -> !(k in ['sample_id', 'run_id', 'fastq_1', 'fastq_2']) }],
                          file(r.fastq_1, checkIfExists: true), file(r.fastq_2, checkIfExists: true)) }

    // 00 input validation; a FAIL stops that sample here
    INPUT_CHECK(reads, thresholds)
    s00 = INPUT_CHECK.out.summary
    passed = s00.filter { meta, js -> status(js) != 'FAIL' }
        .join(reads.map { meta, r1, r2 -> tuple(meta, r1, r2) })
        .map { meta, js, r1, r2 -> tuple(meta, r1, r2) }

    // 01 read QC and trimming
    FASTP(passed, TH['01_read_qc']['fastp'])
    READ_QC_SUMMARY(FASTP.out.report, thresholds)
    s01 = READ_QC_SUMMARY.out.summary
    // a FAIL at 01 stops the sample
    trimmed = s01.filter { meta, js, png -> status(js) != 'FAIL' }
        .join(FASTP.out.reads)
        .map { meta, js, png, r1, r2 -> tuple(meta, r1, r2) }

    // 02 species (runs alongside alignment)
    MYKROBE(trimmed)
    SYLPH(trimmed)
    KRAKEN2(trimmed)
    SPECIES_SUMMARY(MYKROBE.out.json.join(SYLPH.out.tsv).join(KRAKEN2.out.report), thresholds)
    s02 = SPECIES_SUMMARY.out.summary

    // 03 alignment, 04 coverage
    ALIGN(trimmed)
    ALIGN_SUMMARY(ALIGN.out.stats, thresholds)
    s03 = ALIGN_SUMMARY.out.summary
    MOSDEPTH(ALIGN.out.bam)
    COVERAGE_SUMMARY(MOSDEPTH.out.cov, thresholds)
    s04 = COVERAGE_SUMMARY.out.summary

    // 06 deletions: assemble only when a resistance gene drops out (or --always_assemble)
    asm_in = s04.map { meta, js, p1, p2 -> tuple(meta, js) }.join(MOSDEPTH.out.perbase).join(trimmed)
        .branch { meta, js, pb, r1, r2 ->
            assemble: params.always_assemble || dropouts(js)
            skip: true }
    SHOVILL(asm_in.assemble.map { meta, js, pb, r1, r2 -> tuple(meta, r1, r2) })
    CONTIG_ALIGN(SHOVILL.out.contigs)
    DELETION_SUMMARY(
        asm_in.assemble.map { meta, js, pb, r1, r2 -> tuple(meta, js, pb) }
            .join(CONTIG_ALIGN.out.paf.join(SHOVILL.out.contigs))
            .map { meta, js, pb, paf, ctg -> tuple(meta, js, pb, [paf, ctg]) }
            .mix(asm_in.skip.map { meta, js, pb, r1, r2 -> tuple(meta, js, pb, []) }),
        thresholds)
    s06 = DELETION_SUMMARY.out.summary

    // 05 variants
    BCFTOOLS_CALL(ALIGN.out.bam, TH['05_variants'])
    LOFREQ(ALIGN.out.bam, TH['05_variants'])
    FREEBAYES(ALIGN.out.bam, TH['05_variants'])
    NORMALISE(BCFTOOLS_CALL.out.vcf.join(LOFREQ.out.vcf).join(FREEBAYES.out.vcf))
    VARIANT_MERGE(NORMALISE.out.vcfs, thresholds)
    s05 = VARIANT_MERGE.out.summary

    // 07 annotation
    CSQ(VARIANT_MERGE.out.variants)
    SNPEFF(VARIANT_MERGE.out.variants)
    ANNOTATION_SUMMARY(VARIANT_MERGE.out.variants.join(CSQ.out.vcf).join(SNPEFF.out.vcf), thresholds)
    s07 = ANNOTATION_SUMMARY.out.summary

    // 08 lineage
    TBPROFILER(trimmed)
    LINEAGE_SUMMARY(TBPROFILER.out.json, thresholds)
    s08 = LINEAGE_SUMMARY.out.summary

    // 09 resistance grading
    GRADE_RESISTANCE(
        ANNOTATION_SUMMARY.out.annotated
            .join(DELETION_SUMMARY.out.deletions)
            .join(s04.map { meta, js, p1, p2 -> tuple(meta, js) })
            .join(MOSDEPTH.out.perbase)
            .join(TBPROFILER.out.json),
        file("${projectDir}/assets/drugs17.tsv"), thresholds)
    s09 = GRADE_RESISTANCE.out.summary

    // 10 mixed infection, 11 classification
    HET_SITES(ALIGN.out.bam, TH['10_mixed_infection'])
    MIXED_SUMMARY(HET_SITES.out.vcf.join(TBPROFILER.out.json), thresholds)
    s10 = MIXED_SUMMARY.out.summary
    CLASSIFY(GRADE_RESISTANCE.out.results.map { meta, dr, ev -> tuple(meta, dr) }
                 .join(s10.map { meta, js, png -> tuple(meta, js) }), thresholds)
    s11 = CLASSIFY.out.summary

    // step PDFs
    summaries = s00
        .mix(s01.map { meta, js, png -> tuple(meta, js) }, s02, s03,
             s04.map { meta, js, p1, p2 -> tuple(meta, js) }, s05.map { meta, js, png -> tuple(meta, js) },
             s06, s07, s08, s09, s10.map { meta, js, png -> tuple(meta, js) }, s11)

    STEP_REPORT(
        s00.map { meta, js -> tuple(meta, '00_input', js, []) }
           .mix(s01.map { meta, js, png -> tuple(meta, '01_read_qc', js, png) })
           .mix(s02.map { meta, js -> tuple(meta, '02_species', js, []) })
           .mix(s03.map { meta, js -> tuple(meta, '03_alignment', js, []) })
           .mix(s04.map { meta, js, p1, p2 -> tuple(meta, '04_coverage', js, [p1, p2]) })
           .mix(s05.map { meta, js, png -> tuple(meta, '05_variants', js, png) })
           .mix(s06.map { meta, js -> tuple(meta, '06_deletions', js, []) })
           .mix(s07.map { meta, js -> tuple(meta, '07_annotation', js, []) })
           .mix(s08.map { meta, js -> tuple(meta, '08_lineage', js, []) })
           .mix(s09.map { meta, js -> tuple(meta, '09_resistance', js, []) })
           .mix(s10.map { meta, js, png -> tuple(meta, '10_mixed_infection', js, png) })
           .mix(s11.map { meta, js -> tuple(meta, '11_classification', js, []) }),
        step_texts)

    // 12 clinical report (after all 12 step PDFs exist)
    CLINICAL_REPORT(
        summaries.groupTuple(size: 12)
            .join(GRADE_RESISTANCE.out.results)
            .join(ANNOTATION_SUMMARY.out.annotated)
            .join(STEP_REPORT.out.pdf.map { meta, step, pdf -> tuple(meta, pdf) }.groupTuple(size: 12)),
        file("${projectDir}/assets/drugs17.tsv"))

    // R1 run-level QC (all samples of the run)
    run_id = reads.map { meta, r1, r2 -> meta.run_id }.first()
    MULTIQC(run_id,
        FASTP.out.report.map { meta, js, cmd, ver -> js }
            .mix(ALIGN.out.stats.map { meta, fs, md, cmd, ver -> fs })
            .mix(MOSDEPTH.out.cov.map { meta, sm, dist, pb, g, gt, k -> [sm, dist] }.flatten())
            .mix(KRAKEN2.out.report.map { meta, r -> r })
            .collect())
    RUN_SUMMARY(run_id,
        summaries.map { meta, js -> js }.mix(CLINICAL_REPORT.out.summary.map { meta, js -> js }).collect(),
        step_texts)

    // R2 relatedness (on request)
    if (params.relatedness) {
        CONSENSUS(NORMALISE.out.genome.join(MOSDEPTH.out.perbase))
        SNP_DISTS(run_id, CONSENSUS.out.fasta.collect())
        n = CONSENSUS.out.fasta.count()
        IQTREE(run_id, SNP_DISTS.out.out.combine(n).filter { a, m, k -> k >= TH['R2_relatedness']['min_samples_tree'] }
                                       .map { a, m, k -> tuple(a, m) })
        RELATEDNESS_REPORT(run_id, SNP_DISTS.out.out,
            IQTREE.out.ifEmpty(file("${projectDir}/assets/NO_TREE")), thresholds, step_texts)
    }
}
