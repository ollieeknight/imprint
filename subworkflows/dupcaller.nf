include { DUPCALLER_CALL; DUPCALLER_ESTIMATE; DUPCALLER_SUMMARIZE } from '../modules/dupcaller'
include { VEP } from './annotation'
include { fastaRef; indexed; listParam; effectiveMaxZeroQualFraction } from './common'

workflow DUPCALLER {
    take:
        ch_paired_bams // [meta, tumor_bam, tumor_bai, normal_bam, normal_bai]

    main:
        def fasta = fastaRef()
        def h5_indexes = [params.dupcaller_ref_h5, params.dupcaller_tn_h5, params.dupcaller_hp_h5, params.dupcaller_str_h5, params.dupcaller_dbs_h5]
            .collect { h5 -> file(h5) }
        def noise_masks = listParam(params.dupcaller_noise_masks)

        DUPCALLER_CALL(
            ch_paired_bams,
            fasta,
            h5_indexes,
            indexed(params.dupcaller_germline_vcf),
            noise_masks.collect { mask -> file(mask) },
            noise_masks.collect { mask -> file("${mask}.tbi") },
            indexed(params.dupcaller_indel_epon),
            file(params.intervals_bed),
            effectiveMaxZeroQualFraction()
        )
        DUPCALLER_ESTIMATE(DUPCALLER_CALL.out.calls, fasta, h5_indexes)
        DUPCALLER_SUMMARIZE(DUPCALLER_ESTIMATE.out.sample_dir.map { _meta, dir -> dir }.collect(sort: { a, b -> a.name <=> b.name }))

        VEP(
            DUPCALLER_CALL.out.calls.flatMap { meta, dir ->
                ['sbs': "SBS/${meta.pair_id}_sbs.vcf.gz", 'indel': "INDEL/${meta.pair_id}_indel.vcf.gz"].collect { mutation_type, relative ->
                    [meta, mutation_type, dir.resolve(relative), dir.resolve("${relative}.tbi")]
                }
            }
        )

    emit:
        annotated      = VEP.out
        calls          = DUPCALLER_CALL.out.calls
        burden         = DUPCALLER_ESTIMATE.out.burden
        cohort_summary = DUPCALLER_SUMMARIZE.out.summary
        cohort_sbs96   = DUPCALLER_SUMMARIZE.out.sbs96
}
