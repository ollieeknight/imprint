include { MUTECT2_NORMAL_ONLY; GENOMICSDB_IMPORT_PON; CREATE_PON } from '../modules/pon'
include { fastaRef; indexed; callingBed } from './common'

workflow PON_GENERATION {
    take:
        ch_normal_bams // [meta, bam, bai]
        normal_count   // from the samplesheet, so the branch is chosen before any task runs

    main:
        def ch_pon        = channel.empty()
        def ch_normal_vcf = channel.empty()

        if (normal_count >= 2) {
            def fasta    = fastaRef()
            def bed      = callingBed()
            def germline = indexed(params.gnomad_germline_resource_vcf)

            MUTECT2_NORMAL_ONLY(ch_normal_bams, fasta, bed, germline)
            // Sorted so GenomicsDBImport sees the same -V order on every run.
            def ch_vcfs = MUTECT2_NORMAL_ONLY.out.vcf
                .map { _meta, vcf, tbi -> [vcf, tbi] }
                .toSortedList { a, b -> a[0].name <=> b[0].name }
                .map { pairs -> pairs.transpose() }
            GENOMICSDB_IMPORT_PON(ch_vcfs, fasta, bed)
            CREATE_PON(GENOMICSDB_IMPORT_PON.out.genomicsdb, fasta, germline)

            ch_pon        = CREATE_PON.out.pon
            ch_normal_vcf = MUTECT2_NORMAL_ONLY.out.vcf
        } else {
            log.warn "PON_GENERATION: only ${normal_count} normal sample(s); skipping PoN (requires 2 or more)"
        }

    emit:
        pon        = ch_pon
        normal_vcf = ch_normal_vcf
}
