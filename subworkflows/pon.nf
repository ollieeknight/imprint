include { MUTECT2_NORMAL_ONLY  } from '../modules/pon'
include { GENOMICSDB_IMPORT_PON } from '../modules/pon'
include { CREATE_PON           } from '../modules/pon'

workflow PON_GENERATION {
    take:
    ch_normal_bams  // tuple val(meta), path(cram), path(crai)
    normal_count    // int computed at parse time; avoids toList() firing early on pipeline error

    main:
    ch_pon = channel.empty()
    ch_normal_vcf = channel.empty()

    if (normal_count >= 2) {
        MUTECT2_NORMAL_ONLY(ch_normal_bams)

        def ch_normal_vcf_split = MUTECT2_NORMAL_ONLY.out.vcf.multiMap { meta, vcf, tbi ->
            pon_input: [meta, vcf, tbi]
            manifest:  [meta, vcf, tbi]
        }

        def ch_pon_input = ch_normal_vcf_split.pon_input
            .map { _meta, vcf, tbi -> ['cohort', vcf, tbi] }
            .groupTuple()

        GENOMICSDB_IMPORT_PON(ch_pon_input)
        CREATE_PON(GENOMICSDB_IMPORT_PON.out.genomicsdb)
        ch_pon = CREATE_PON.out.pon
        ch_normal_vcf = ch_normal_vcf_split.manifest
    } else if (!params.dupcaller) {
        log.warn "PON_GENERATION: only ${normal_count} normal sample(s); skipping PoN (requires ≥2)"
    }

    emit:
    pon        = ch_pon
    normal_vcf = ch_normal_vcf
}
