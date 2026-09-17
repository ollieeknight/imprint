include {
    STRELKA2_MERGE
    VEP_ANNOTATE
    ENSEMBLE_CONSENSUS
    VAFATOR
    TUMOR_NORMAL_FILTER
} from '../modules/annotation'

include { VARLOCIRAPTOR_MERGE } from '../modules/varlociraptor'
include { VARLOCIRAPTOR } from './varlociraptor'

workflow ANNOTATION {
    take:
        ch_paired_bams
        ch_mutect2
        ch_strelka_snv
        ch_strelka_indel
        ch_deepsomatic
        ch_sex

    main:
        ch_strelka_snv
            .map { meta, vcf, tbi -> [meta.pair_id, meta, vcf, tbi] }
            .join(
                ch_strelka_indel.map { meta, vcf, tbi -> [meta.pair_id, vcf, tbi] },
                by: 0, failOnDuplicate: false, failOnMismatch: true
            )
            .map { _pair_id, meta, snv, snv_tbi, indel, indel_tbi -> [meta, snv, snv_tbi, indel, indel_tbi] }
            .set { ch_strelka_for_concat }

        STRELKA2_MERGE(ch_strelka_for_concat)

        ch_mutect2.map { meta, vcf, tbi -> [meta.pair_id, meta, [m2_vcf: vcf, m2_tbi: tbi]] }
            .mix(STRELKA2_MERGE.out.all_vcf.map { meta, vcf, tbi -> [meta.pair_id, meta, [st_vcf: vcf, st_tbi: tbi]] })
            .mix(ch_deepsomatic.map { meta, vcf, tbi -> [meta.pair_id, meta, [ds_vcf: vcf, ds_tbi: tbi]] })
            .groupTuple(by: 0, size: 3)
            .map { pair_id, metas, results ->
                def meta     = metas[0]
                def combined = results.inject([:]) { acc, item -> acc + item }
                assert combined.containsKey('m2_vcf') && combined.containsKey('st_vcf') && combined.containsKey('ds_vcf') :
                    "Ensemble for ${pair_id}: missing mutect2, strelka2 or deepsomatic output"
                [meta, combined.m2_vcf, combined.m2_tbi, combined.st_vcf, combined.st_tbi, combined.ds_vcf, combined.ds_tbi]
            }
            .set { ch_ensemble_input }

        ENSEMBLE_CONSENSUS(ch_ensemble_input)

        def ch_ensemble_split = ENSEMBLE_CONSENSUS.out.vcf.multiMap { meta, vcf, tbi ->
            vafator:       [meta, vcf, tbi]
            varlociraptor: [meta, vcf, tbi]
        }

        ch_ensemble_split.vafator
            .map { meta, vcf, tbi -> [meta.pair_id, meta, vcf, tbi] }
            .join(
                ch_paired_bams.map { meta, tb, tbai, nb, nbai -> [meta.pair_id, meta, tb, tbai, nb, nbai] },
                by: 0, failOnMismatch: true
            )
            .map { _pair_id, meta, vcf, tbi, _paired_meta, tb, tbai, nb, nbai ->
                [meta, vcf, tbi, tb, tbai, nb, nbai]
            }
            .set { ch_vafator_input }

        VAFATOR(ch_vafator_input)

        VARLOCIRAPTOR(ch_paired_bams, ch_ensemble_split.varlociraptor, ch_sex)

        TUMOR_NORMAL_FILTER(VAFATOR.out.vcf)

        VEP_ANNOTATE(
            TUMOR_NORMAL_FILTER.out.vcf.map { meta, vcf, tbi -> [meta, 'somatic', vcf, tbi] },
            file(params.spliceai_snv_vcf),
            file("${params.spliceai_snv_vcf}.tbi"),
            file(params.spliceai_indel_vcf ?: params.spliceai_snv_vcf),
            file("${params.spliceai_indel_vcf ?: params.spliceai_snv_vcf}.tbi"),
            file(params.dbnsfp_gz),
            file("${params.dbnsfp_gz}.tbi")
        )

        VEP_ANNOTATE.out.vcf
            .map { meta, _mutation_type, vcf, tbi -> [meta.pair_id, meta, vcf, tbi] }
            .join(
                VARLOCIRAPTOR.out.map { meta, bcf, csi -> [meta.pair_id, bcf, csi] },
                by: 0, failOnMismatch: true
            )
            .map { _pair_id, meta, vcf, tbi, bcf, csi -> [meta, vcf, tbi, bcf, csi] }
            .set { ch_varlociraptor_merge_input }

        VARLOCIRAPTOR_MERGE(ch_varlociraptor_merge_input)

    emit:
        final_vcf     = VARLOCIRAPTOR_MERGE.out.vcf
        varlociraptor = VARLOCIRAPTOR.out // [meta, bcf, csi]
}
