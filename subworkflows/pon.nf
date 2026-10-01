include { MUTECT2_NORMAL_ONLY; GENOMICSDB_IMPORT_PON; CREATE_PON } from '../modules/pon'

workflow PON_GENERATION {
    take:
    ch_normal_bams  // tuple val(meta), path(cram), path(crai)
    normal_count    // from the samplesheet, so the PoN branch is chosen before any task runs

    main:
    ch_pon = channel.empty()
    ch_normal_vcf = channel.empty()

    if (normal_count >= 2) {
        MUTECT2_NORMAL_ONLY(ch_normal_bams)

        // Sorted so GenomicsDBImport sees the same -V order on every run.
        def ch_pon_input = MUTECT2_NORMAL_ONLY.out.vcf
            .map { _meta, vcf, tbi -> ['cohort', vcf, tbi] }
            .groupTuple()
            .map { cohort, vcfs, tbis -> [cohort, vcfs.sort { v -> v.name }, tbis.sort { t -> t.name }] }

        GENOMICSDB_IMPORT_PON(ch_pon_input)
        CREATE_PON(GENOMICSDB_IMPORT_PON.out.genomicsdb)
        ch_pon = CREATE_PON.out.pon
        ch_normal_vcf = MUTECT2_NORMAL_ONLY.out.vcf
    } else if (!params.dupcaller) {
        log.warn "PON_GENERATION: only ${normal_count} normal sample(s); skipping PoN (requires ≥2)"
    }

    emit:
    pon        = ch_pon
    normal_vcf = ch_normal_vcf
}
