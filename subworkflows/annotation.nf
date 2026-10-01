include { STRELKA2_MERGE; ENSEMBLE_CONSENSUS; VAFATOR; MARK_ON_TARGET; VEP_ANNOTATE } from '../modules/annotation'
include { VARLOCIRAPTOR_MERGE } from '../modules/varlociraptor'
include { VARLOCIRAPTOR }       from './varlociraptor'
include { fastaRef; indexed }   from './common'

// VEP with its references; DUPCALLER uses it too.
workflow VEP {
    take:
        ch_vcfs // [meta, mutation_type, vcf, tbi]

    main:
        def spliceai_indel = params.spliceai_indel_vcf ?: params.spliceai_snv_vcf
        VEP_ANNOTATE(
            ch_vcfs,
            fastaRef(),
            [file(params.vep_cache), file(params.vep_plugins_dir)],
            indexed(params.spliceai_snv_vcf) + indexed(spliceai_indel),
            indexed(params.dbnsfp_gz),
            indexed(params.alphamissense_tsv),
            indexed(params.gnomad_exomes_vep_vcf) + indexed(params.gnomad_genomes_vep_vcf),
            indexed(params.cosmic_vcf) + indexed(params.cosmic_noncoding_vcf)
        )

    emit:
        VEP_ANNOTATE.out.vcf
}

workflow ANNOTATION {
    take:
        ch_paired_bams // [meta, tumor_bam, tumor_bai, normal_bam, normal_bai], the calling BAMs
        ch_mutect2
        ch_strelka_snv
        ch_strelka_indel
        ch_deepsomatic
        ch_sex

    main:
        STRELKA2_MERGE(
            ch_strelka_snv
                .map { meta, vcf, tbi -> [meta.pair_id, meta, vcf, tbi] }
                .join(ch_strelka_indel.map { meta, vcf, tbi -> [meta.pair_id, vcf, tbi] }, failOnMismatch: true)
                .map { _pair_id, meta, snv, snv_tbi, indel, indel_tbi -> [meta, snv, snv_tbi, indel, indel_tbi] }
        )

        def ch_ensemble_input = ch_mutect2
            .map { meta, vcf, tbi -> [meta.pair_id, meta, vcf, tbi] }
            .join(STRELKA2_MERGE.out.all_vcf.map { meta, vcf, tbi -> [meta.pair_id, vcf, tbi] }, failOnMismatch: true)
            .join(ch_deepsomatic.map { meta, vcf, tbi -> [meta.pair_id, vcf, tbi] }, failOnMismatch: true)
            .map { _pair_id, meta, m2, m2_tbi, st, st_tbi, ds, ds_tbi -> [meta, m2, m2_tbi, st, st_tbi, ds, ds_tbi] }
        ENSEMBLE_CONSENSUS(ch_ensemble_input, fastaRef())

        VAFATOR(
            ENSEMBLE_CONSENSUS.out.vcf
                .map { meta, vcf, tbi -> [meta.pair_id, meta, vcf, tbi] }
                .join(ch_paired_bams.map { meta, tb, tbai, nb, nbai -> [meta.pair_id, tb, tbai, nb, nbai] }, failOnMismatch: true)
                .map { _pair_id, meta, vcf, tbi, tb, tbai, nb, nbai -> [meta, vcf, tbi, tb, tbai, nb, nbai] }
        )

        VARLOCIRAPTOR(ch_paired_bams, ENSEMBLE_CONSENSUS.out.vcf, ch_sex)

        // Off-target runs call on padded targets; flag sites inside the unpadded ones.
        MARK_ON_TARGET(VAFATOR.out.vcf, params.off_target ? file(params.intervals_bed) : [])

        VEP(MARK_ON_TARGET.out.vcf.map { meta, vcf, tbi -> [meta, 'somatic', vcf, tbi] })

        VARLOCIRAPTOR_MERGE(
            VEP.out
                .map { meta, _mutation_type, vcf, tbi -> [meta.pair_id, meta, vcf, tbi] }
                .join(VARLOCIRAPTOR.out.map { meta, bcf, csi -> [meta.pair_id, bcf, csi] }, failOnMismatch: true)
                .map { _pair_id, meta, vcf, tbi, bcf, csi -> [meta, vcf, tbi, bcf, csi] }
        )

    emit:
        final_vcf     = VARLOCIRAPTOR_MERGE.out.vcf
        varlociraptor = VARLOCIRAPTOR.out
}
