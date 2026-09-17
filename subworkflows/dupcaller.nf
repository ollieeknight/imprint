include { DUPCALLER_CALL; DUPCALLER_ESTIMATE; DUPCALLER_SUMMARIZE } from '../modules/dupcaller'
include { VEP_ANNOTATE } from '../modules/annotation'

workflow DUPCALLER {
    take:
        ch_paired_bams

    main:
        if (params.dupcaller) {
            def noiseMasks = params.dupcaller_noise_masks instanceof Collection
                ? params.dupcaller_noise_masks.findAll { item -> item }
                : (params.dupcaller_noise_masks ? params.dupcaller_noise_masks.toString().tokenize(',').collect { item -> item.trim() }.findAll { item -> item } : [])
            def stagedNoiseMasks = noiseMasks ? noiseMasks.collect { mask -> file(mask) } : [file(params.dupcaller_umi_allowlist)]
            def stagedNoiseIndexes = noiseMasks ? noiseMasks.collect { mask -> file("${mask}.tbi") } : [file("${projectDir}/bin/validate_dupcaller_tags.awk")]
            def maxZeroQualFraction = params.dupcaller_max_zero_qual_fraction != null
                ? params.dupcaller_max_zero_qual_fraction
                : (noiseMasks ? 0.5 : 0.1)

            DUPCALLER_CALL(
                ch_paired_bams,
                file(params.ref_fasta),
                file(params.genome_fai),
                file(params.dupcaller_ref_h5),
                file(params.dupcaller_tn_h5),
                file(params.dupcaller_hp_h5),
                file(params.dupcaller_str_h5),
                file(params.dupcaller_dbs_h5),
                file(params.dupcaller_germline_vcf),
                file("${params.dupcaller_germline_vcf}.tbi"),
                stagedNoiseMasks,
                stagedNoiseIndexes,
                params.dupcaller_indel_epon ? file(params.dupcaller_indel_epon) : file(params.dupcaller_umi_allowlist),
                params.dupcaller_indel_epon ? file("${params.dupcaller_indel_epon}.tbi") : file("${projectDir}/bin/validate_dupcaller_tags.awk"),
                file(params.intervals_bed),
                noiseMasks as boolean,
                maxZeroQualFraction
            )

            DUPCALLER_ESTIMATE(
                DUPCALLER_CALL.out.calls,
                file(params.ref_fasta),
                file(params.dupcaller_ref_h5),
                file(params.dupcaller_tn_h5),
                file(params.dupcaller_hp_h5),
                file(params.dupcaller_str_h5),
                file(params.dupcaller_dbs_h5)
            )

            DUPCALLER_SUMMARIZE(
                DUPCALLER_ESTIMATE.out.sample_dir.map { _meta, dir -> dir }.collect()
            )

            def ch_vep_input = DUPCALLER_CALL.out.calls.flatMap { meta, dir ->
                [['sbs', "SBS/${meta.pair_id}_sbs.vcf.gz"],
                 ['indel', "INDEL/${meta.pair_id}_indel.vcf.gz"]].collect { mutation_type, relative ->
                    [meta, mutation_type, dir.resolve(relative), dir.resolve("${relative}.tbi")]
                }
            }

            VEP_ANNOTATE(
                ch_vep_input,
                file(params.spliceai_snv_vcf),
                file("${params.spliceai_snv_vcf}.tbi"),
                file(params.spliceai_indel_vcf ?: params.spliceai_snv_vcf),
                file("${params.spliceai_indel_vcf ?: params.spliceai_snv_vcf}.tbi"),
                file(params.dbnsfp_gz),
                file("${params.dbnsfp_gz}.tbi")
            )

            ch_annotated = VEP_ANNOTATE.out.vcf
            ch_calls = DUPCALLER_CALL.out.calls
            ch_burden = DUPCALLER_ESTIMATE.out.burden
            ch_summary = DUPCALLER_SUMMARIZE.out.summary
            ch_sbs96 = DUPCALLER_SUMMARIZE.out.sbs96
        } else {
            ch_annotated = channel.empty()
            ch_calls = channel.empty()
            ch_burden = channel.empty()
            ch_summary = channel.empty()
            ch_sbs96 = channel.empty()
        }

    emit:
        annotated = ch_annotated
        calls = ch_calls
        burden = ch_burden
        cohort_summary = ch_summary
        cohort_sbs96 = ch_sbs96
}
