include {
    VARLOCIRAPTOR_ALIGNMENT_PROPERTIES
    VARLOCIRAPTOR_PREPROCESS
    VARLOCIRAPTOR_CALLVARIANTS
    VARLOCIRAPTOR_INDEX
} from '../modules/varlociraptor'
include { fastaRef; tumors; normalsOnce } from './common'

workflow VARLOCIRAPTOR {
    take:
        ch_paired_bams // [meta, tumor_bam, tumor_bai, normal_bam, normal_bai]
        ch_candidates  // [meta, vcf, tbi], the consensus VCF
        ch_sex         // [[sample_id, sex], ...] as one list

    main:
        def fasta = fastaRef()

        // Somalier's sex call for the normal picks the scenario; UNKNOWN uses xx.
        def ch_scenario = ch_paired_bams
            .combine(ch_sex.map { sexes -> [sexes] })
            .map { meta, _tb, _tbai, _nb, _nbai, sexes ->
                def xy = sexes.collectEntries { s -> s }[meta.normal_id] == 'XY'
                [meta.pair_id, file(xy ? params.varlociraptor_scenario_xy : params.varlociraptor_scenario_xx)]
            }

        // Alignment properties once per sample; preprocessing per pair, since each pair has its own candidates.
        VARLOCIRAPTOR_ALIGNMENT_PROPERTIES(tumors(ch_paired_bams).mix(normalsOnce(ch_paired_bams)), fasta)

        def ch_preprocess_input = tumors(ch_paired_bams)
            .mix(ch_paired_bams.map { meta, _tb, _tbai, nb, nbai -> [meta, 'normal', meta.normal_id, nb, nbai] })
            .map { meta, role, sample_id, bam, bai -> [meta.pair_id, meta, role, sample_id, bam, bai] }
            .combine(ch_candidates.map { meta, vcf, tbi -> [meta.pair_id, vcf, tbi] }, by: 0)
            .map { _pair_id, meta, role, sample_id, bam, bai, vcf, tbi -> [sample_id, meta, role, bam, bai, vcf, tbi] }
            .combine(VARLOCIRAPTOR_ALIGNMENT_PROPERTIES.out.json.map { _meta, _role, sample_id, json -> [sample_id, json] }, by: 0)
            .map { sample_id, meta, role, bam, bai, vcf, tbi, json -> [meta, role, sample_id, bam, bai, vcf, tbi, json] }
        VARLOCIRAPTOR_PREPROCESS(ch_preprocess_input, fasta)

        def ch_obs = VARLOCIRAPTOR_PREPROCESS.out.bcf.branch { _meta, role, _bcf ->
            tumor: role == 'tumor'
            normal: role == 'normal'
        }
        def ch_call_input = ch_obs.tumor
            .map { meta, _role, bcf -> [meta.pair_id, meta, bcf] }
            .join(ch_obs.normal.map { meta, _role, bcf -> [meta.pair_id, bcf] }, failOnMismatch: true)
            .join(ch_scenario, failOnMismatch: true)
            .map { _pair_id, meta, tumor_bcf, normal_bcf, scenario -> [meta, tumor_bcf, normal_bcf, scenario] }
        VARLOCIRAPTOR_CALLVARIANTS(ch_call_input)
        VARLOCIRAPTOR_INDEX(VARLOCIRAPTOR_CALLVARIANTS.out.bcf)

    emit:
        VARLOCIRAPTOR_INDEX.out.bcf // [meta, bcf, csi]
}
