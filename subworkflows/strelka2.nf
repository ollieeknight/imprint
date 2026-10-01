include { MANTA_PREP_BAM; MANTA; STRELKA2_CALL } from '../modules/strelka2'
include { fastaRef; tumors; normalsOnce; repair } from './common'

workflow STRELKA2 {
    take:
        ch_paired_bams  // [meta, tumor_bam, tumor_bai, normal_bam, normal_bai]
        ch_intervals_gz // [intervals_gz, intervals_tbi], or [[], []] for WGS

    main:
        def fasta = fastaRef()

        MANTA_PREP_BAM(tumors(ch_paired_bams).mix(normalsOnce(ch_paired_bams)))
        def ch_prepped = repair(MANTA_PREP_BAM.out.bam)

        MANTA(ch_prepped, ch_intervals_gz, fasta)

        def ch_strelka_input = ch_prepped
            .map { meta, tb, tbai, nb, nbai -> [meta.pair_id, meta, tb, tbai, nb, nbai] }
            .join(MANTA.out.indels.map { meta, indels, tbi -> [meta.pair_id, indels, tbi] }, failOnDuplicate: true, failOnMismatch: true)
            .map { _pair_id, meta, tb, tbai, nb, nbai, indels, indels_tbi -> [meta, tb, tbai, nb, nbai, indels, indels_tbi] }
        STRELKA2_CALL(ch_strelka_input, ch_intervals_gz, fasta)

    emit:
        snv      = STRELKA2_CALL.out.vcf
        indel    = STRELKA2_CALL.out.indels
        manta_sv = MANTA.out.svs
}
