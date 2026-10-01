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

        def ch_tumours = ch_paired_bams.map { meta, tb, tbai, _nb, _nbai -> [meta, 'tumour', meta.tumor_id, tb, tbai] }

        // Alignment properties run once per normal, not per pair; preprocessing stays
        // per pair because each pair has its own candidates.
        def ch_normals = ch_paired_bams
            .map { meta, _tb, _tbai, nb, nbai -> [groupKey(meta.normal_id, meta.tumour_count), meta, nb, nbai] }
            .groupTuple(by: 0)
            .map { normal_id, metas, nbs, nbais -> [metas.min { m -> m.pair_id }, 'normal', normal_id.toString(), nbs[0], nbais[0]] }

        VARLOCIRAPTOR_ALIGNMENT_PROPERTIES(ch_tumours.mix(ch_normals))

        def ch_preprocess_input = ch_tumours
            .mix(ch_paired_bams.map { meta, _tb, _tbai, nb, nbai -> [meta, 'normal', meta.normal_id, nb, nbai] })
            .map { meta, role, sample_id, bam, bai -> [meta.pair_id, meta, role, sample_id, bam, bai] }
            .combine(ch_candidates.map { meta, vcf, tbi -> [meta.pair_id, vcf, tbi] }, by: 0)
            .map { _pair_id, meta, role, sample_id, bam, bai, vcf, tbi -> [sample_id, meta, role, bam, bai, vcf, tbi] }
            .combine(VARLOCIRAPTOR_ALIGNMENT_PROPERTIES.out.json.map { _meta, _role, sample_id, json -> [sample_id, json] }, by: 0)
            .map { sample_id, meta, role, bam, bai, vcf, tbi, json -> [meta, role, sample_id, bam, bai, vcf, tbi, json] }

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
