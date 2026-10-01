include { PREPARE_INTERVALS; SPLIT_INTERVALS } from '../modules/intervals'
include { BASE_RECALIBRATOR; APPLY_BQSR }      from '../modules/bqsr'
include { DEEPSOMATIC }                        from '../modules/deepsomatic'
include { MUTECT2 }                            from './mutect2'
include { STRELKA2 }                           from './strelka2'
include { fastaRef; indexed; callingBed; tumors; normalsOnce; repair } from './common'

workflow SOMATIC_CALLING {
    take:
        ch_paired_bams // [meta, tumor_bam, tumor_bai, normal_bam, normal_bai]

    main:
        def fasta = fastaRef()
        def bed   = callingBed()

        def ch_intervals_gz = channel.value([[], []])
        if (bed) {
            PREPARE_INTERVALS(bed)
            ch_intervals_gz = PREPARE_INTERVALS.out.intervals
        }

        SPLIT_INTERVALS(bed, fasta)
        def ch_shards = SPLIT_INTERVALS.out.intervals_gz
            .flatten()
            .map { f -> [f.simpleName, f] }
            .join(SPLIT_INTERVALS.out.intervals_tbi.flatten().map { f -> [f.simpleName, f] }, failOnDuplicate: true, failOnMismatch: true)
            .map { _name, gz, tbi -> [gz, tbi] }

        def ch_calling_bams = ch_paired_bams
        def ch_bqsr_recal   = channel.empty()
        if (!params.skip_bqsr) {
            def ch_samples = tumors(ch_paired_bams).mix(normalsOnce(ch_paired_bams))
            def known_sites = indexed(params.dbsnp) + indexed(params.known_indels_mills) + indexed(params.known_snps_1000g)
            BASE_RECALIBRATOR(ch_samples, fasta, bed, known_sites)

            def ch_apply_input = ch_samples
                .map { meta, role, sample_id, bam, bai -> [sample_id, meta, role, bam, bai] }
                .join(BASE_RECALIBRATOR.out.recal_table.map { _meta, _role, sample_id, table -> [sample_id, table] }, failOnDuplicate: true, failOnMismatch: true)
                .map { sample_id, meta, role, bam, bai, table -> [meta, role, sample_id, bam, bai, table] }
            APPLY_BQSR(ch_apply_input, fasta)

            ch_calling_bams = repair(APPLY_BQSR.out.bam)
            ch_bqsr_recal   = BASE_RECALIBRATOR.out.recal_table
        }

        MUTECT2(ch_calling_bams, ch_shards)
        STRELKA2(ch_calling_bams, ch_intervals_gz)
        DEEPSOMATIC(ch_calling_bams, fasta, bed)

    emit:
        calling_bams          = ch_calling_bams
        mutect2_vcf           = MUTECT2.out.vcf
        mutect2_raw_vcf       = MUTECT2.out.raw_vcf
        mutect2_contamination = MUTECT2.out.contamination
        mutect2_segments      = MUTECT2.out.segments
        strelka_snv           = STRELKA2.out.snv
        strelka_indel         = STRELKA2.out.indel
        manta_sv              = STRELKA2.out.manta_sv
        deepsomatic_vcf       = DEEPSOMATIC.out.vcf
        bqsr_recal            = ch_bqsr_recal
}
