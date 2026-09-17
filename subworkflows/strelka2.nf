include { PREP_MANTA_BAM }    from '../modules/align'
include { MANTA; CALL }              from '../modules/strelka2'

workflow STRELKA2 {
    take:
        ch_paired_bams  // [meta, tumor_bam, tumor_bai, normal_bam, normal_bai]
        ch_intervals_gz // [intervals_gz, intervals_tbi]

    main:
        ch_paired_bams.multiMap { meta, tb, tbai, nb, nbai ->
            tumour: [meta, 'tumour', meta.tumor_id,  tb, tbai]
            normal: [meta, 'normal', meta.normal_id, nb, nbai]
        }.set { ch_split }

        PREP_MANTA_BAM(ch_split.tumour.mix(ch_split.normal))

        PREP_MANTA_BAM.out.bam
            .branch { _meta, role, _sample_id, _bam, _bai ->
                tumour: role == 'tumour'
                normal: role == 'normal'
            }
            .set { ch_stripped }

        ch_stripped.tumour
            .map { meta, _role, _sample_id, tb, tbai -> [meta.pair_id, meta, tb, tbai] }
            .join(
                ch_stripped.normal.map { meta, _role, _sample_id, nb, nbai -> [meta.pair_id, nb, nbai] },
                by: 0, failOnDuplicate: true, failOnMismatch: true
            )
            .map { _pair_id, meta, tb, tbai, nb, nbai -> [meta, tb, tbai, nb, nbai] }
            .set { ch_stripped_paired }

        MANTA(ch_stripped_paired, ch_intervals_gz)

        ch_stripped_paired
            .map { meta, tb, tbai, nb, nbai -> [meta.pair_id, meta, tb, tbai, nb, nbai] }
            .join(
                MANTA.out.indels.map { meta, indels, tbi -> [meta.pair_id, indels, tbi] },
                by: 0, failOnDuplicate: true, failOnMismatch: true
            )
            .map { _pair_id, meta, tb, tbai, nb, nbai, indels, indels_tbi ->
                [meta, tb, tbai, nb, nbai, indels, indels_tbi]
            }
            .set { ch_strelka_input }

        CALL(ch_strelka_input, ch_intervals_gz)

    emit:
        snv   = CALL.out.vcf
        indel = CALL.out.indels
        manta_sv = MANTA.out.svs
}
