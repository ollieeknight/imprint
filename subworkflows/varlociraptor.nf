include {
    VARLOCIRAPTOR_ALIGNMENT_PROPERTIES
    VARLOCIRAPTOR_PREPROCESS
    VARLOCIRAPTOR_CALLVARIANTS
    VARLOCIRAPTOR_INDEX
} from '../modules/varlociraptor'

workflow VARLOCIRAPTOR {
    take:
        ch_paired_bams // [meta, tumor_bam, tumor_bai, normal_bam, normal_bai]
        ch_candidates  // [meta, vcf, tbi], pre-VAFATOR/VEP consensus VCF
        ch_sex         // [[sample_id, sex], ...] single list (same channel CHARACTERISATION uses)

    main:
        // Select the ploidy scenario from Somalier sex; UNKNOWN uses XX.
        def ch_scenario = ch_paired_bams
            .combine(ch_sex.map { it -> [it] })
            .map { meta, _tb, _tbai, _nb, _nbai, sex_list ->
                def sex_map = sex_list.collectEntries { v -> v }
                def sex = sex_map[meta.normal_id]
                def scenario = (sex == 'XY')
                    ? file(params.varlociraptor_scenario_xy)
                    : file(params.varlociraptor_scenario_xx)
                [meta.pair_id, scenario]
            }

        // Keep each tumour and normal stream keyed to its pair.
        def ch_sample_bams = ch_paired_bams.flatMap { meta, tb, tbai, nb, nbai ->
            [
                [[meta.pair_id, 'tumour'], meta, 'tumour', meta.tumor_id,  tb, tbai],
                [[meta.pair_id, 'normal'], meta, 'normal', meta.normal_id, nb, nbai],
            ]
        }

        VARLOCIRAPTOR_ALIGNMENT_PROPERTIES(
            ch_sample_bams.map { _key, meta, role, sample_id, bam, bai -> [meta, role, sample_id, bam, bai] }
        )

        def ch_alignment_properties = VARLOCIRAPTOR_ALIGNMENT_PROPERTIES.out.json
            .map { meta, role, _sample_id, json -> [[meta.pair_id, role], json] }

        def ch_preprocess_input = ch_sample_bams
            .map { key, meta, role, sample_id, bam, bai -> [meta.pair_id, key, meta, role, sample_id, bam, bai] }
            .combine(ch_candidates.map { meta, vcf, tbi -> [meta.pair_id, vcf, tbi] }, by: 0)
            .map { _pair_id, key, meta, role, sample_id, bam, bai, vcf, tbi -> [key, meta, role, sample_id, bam, bai, vcf, tbi] }
            .join(ch_alignment_properties, by: 0, failOnMismatch: true)
            .map { _key, meta, role, sample_id, bam, bai, vcf, tbi, json -> [meta, role, sample_id, bam, bai, vcf, tbi, json] }

        VARLOCIRAPTOR_PREPROCESS(ch_preprocess_input)

        def ch_tumour_obs = VARLOCIRAPTOR_PREPROCESS.out.bcf
            .filter { _meta, role, _bcf -> role == 'tumour' }
            .map { meta, _role, bcf -> [meta.pair_id, meta, bcf] }
        def ch_normal_obs = VARLOCIRAPTOR_PREPROCESS.out.bcf
            .filter { _meta, role, _bcf -> role == 'normal' }
            .map { meta, _role, bcf -> [meta.pair_id, bcf] }

        def ch_call_input = ch_tumour_obs
            .join(ch_normal_obs, by: 0, failOnMismatch: true)
            .join(ch_scenario, by: 0, failOnMismatch: true)
            .map { _pair_id, meta, tumour_bcf, normal_bcf, scenario -> [meta, tumour_bcf, normal_bcf, scenario] }

        VARLOCIRAPTOR_CALLVARIANTS(ch_call_input)
        VARLOCIRAPTOR_INDEX(VARLOCIRAPTOR_CALLVARIANTS.out.bcf)

    emit:
        VARLOCIRAPTOR_INDEX.out.bcf // [meta, bcf, csi]
}
