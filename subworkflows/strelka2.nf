include { PREP_MANTA_BAM }    from '../modules/align'
include { MANTA; CALL }              from '../modules/strelka2'

workflow STRELKA2 {
    take:
        ch_paired_bams  // [meta, tumor_bam, tumor_bai, normal_bam, normal_bai]
        ch_intervals_gz // [intervals_gz, intervals_tbi]

    main:
        // One task per normal, not per pair, carrying the first pair's meta so the hash is stable.
        def ch_normals = ch_paired_bams
            .map { meta, _tb, _tbai, nb, nbai -> [groupKey(meta.normal_id, meta.tumour_count), meta, nb, nbai] }
            .groupTuple(by: 0)
            .map { normal_id, metas, nbs, nbais -> [metas.min { m -> m.pair_id }, 'normal', normal_id.toString(), nbs[0], nbais[0]] }

        PREP_MANTA_BAM(
            ch_paired_bams
                .map { meta, tb, tbai, _nb, _nbai -> [meta, 'tumour', meta.tumor_id, tb, tbai] }
                .mix(ch_normals)
        )

        PREP_MANTA_BAM.out.bam
            .branch { _meta, role, _sample_id, _bam, _bai ->
                tumour: role == 'tumour'
                normal: role == 'normal'
            }
            .set { ch_stripped }

        ch_stripped.tumour
            .map { meta, _role, _sample_id, tb, tbai -> [meta.normal_id, meta, tb, tbai] }
            .combine(ch_stripped.normal.map { _meta, _role, normal_id, nb, nbai -> [normal_id, nb, nbai] }, by: 0)
            .map { _normal_id, meta, tb, tbai, nb, nbai -> [meta, tb, tbai, nb, nbai] }
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
