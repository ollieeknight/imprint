include {
    MERGE_DONOR_BAMS
    OPTITYPE_EXTRACT
    OPTITYPE_GENOTYPE
    KIR_MAPPER
    KIR_COLLATE
    PATHSEQ
    TELSEQ
} from '../modules/tumour_characterisation'

include { MIXCR; MIXCR_EXPORT_CLONES } from '../modules/mixcr'

def isExtraEnabled(String name, Set aliases = []) {
    if (params.extras) {
        def _extras = params.extras.tokenize(',').collect { v -> v.trim().toLowerCase() }.toSet()
        return (name in _extras) || aliases.any { a -> a in _extras }
    }
    return false
}

workflow CHARACTERISATION {
    take:
        ch_all_bams_per_donor // [meta_donor, bams, bais], all analysis BAMs per donor
        ch_bams               // [meta, bam, bai], all sample analysis BAMs (for PathSeq)
        ch_trimmed_reads      // [meta, [r1s], [r2s]], trimmed FASTQs (for MiXCR)
        ch_merged_bam         // [meta, bam, bai], pre-dedup merged BAMs (for TelSeq)
        ch_fastp_stats        // [id, short_inserts, read_length], for TelSeq read length
        ch_pathseq_references // [host IMG, host HSS, microbe IMG, microbe dict, taxonomy DB]

    main:
        def run_hla     = isExtraEnabled('hla', ['optitype'] as Set)
        def run_kir     = isExtraEnabled('kir', ['kirmapper', 'kir_mapper'] as Set)
        def run_pathseq = isExtraEnabled('pathseq')
        def run_mixcr   = isExtraEnabled('mixcr')
        def run_telseq  = isExtraEnabled('telseq')

        if (run_hla || run_kir) {
            MERGE_DONOR_BAMS(ch_all_bams_per_donor)
            ch_donor_bams_split = MERGE_DONOR_BAMS.out.bam.multiMap { meta, bam, bai ->
                hla: [meta, bam, bai]
                kir: [meta, bam, bai]
            }
        }

        if (run_hla) {
            OPTITYPE_EXTRACT(ch_donor_bams_split.hla, channel.fromPath("${params.hla_reference}*").collect())
            OPTITYPE_GENOTYPE(OPTITYPE_EXTRACT.out.bams)
        }

        if (run_kir) {
            KIR_MAPPER(ch_donor_bams_split.kir)
            KIR_COLLATE(
                KIR_MAPPER.out.raw,
                file("${projectDir}/bin/collate_kir.py")
            )
        }

        if (run_pathseq) {
            PATHSEQ(ch_bams, ch_pathseq_references)
        }

        if (run_mixcr) {
            def ch_mixcr_license = params.mixcr_license ? file(params.mixcr_license) : file('NO_FILE')
            MIXCR(ch_trimmed_reads, ch_mixcr_license)
            MIXCR_EXPORT_CLONES(MIXCR.out.clns, ch_mixcr_license)
        }


        if (run_telseq) {
            ch_merged_bam
                .map { meta, bam, bai -> [meta.id, meta, bam, bai] }
                .join(ch_fastp_stats)
                .map { _id, meta, bam, bai, _short_inserts, read_length ->
                    [meta + [read_length: read_length], bam, bai]
                }
                .set { ch_bams_with_rlen }
            TELSEQ(ch_bams_with_rlen)
        }

    emit:
        hla_result     = run_hla     ? OPTITYPE_GENOTYPE.out.result : channel.empty()
        hla_plot       = run_hla     ? OPTITYPE_GENOTYPE.out.plot   : channel.empty()
        kir_ncopy      = run_kir     ? KIR_COLLATE.out.copy_number  : channel.empty()
        kir_calls      = run_kir     ? KIR_COLLATE.out.calls        : channel.empty()
        kir_reports    = run_kir     ? KIR_COLLATE.out.candidates   : channel.empty()
        kir_raw_archive = run_kir    ? KIR_COLLATE.out.raw_archive  : channel.empty()
        pathseq_bam    = run_pathseq ? PATHSEQ.out.bam              : channel.empty()
        pathseq_scores = run_pathseq ? PATHSEQ.out.scores           : channel.empty()
        pathseq_filter_metrics = run_pathseq ? PATHSEQ.out.filter_metrics : channel.empty()
        pathseq_score_metrics  = run_pathseq ? PATHSEQ.out.score_metrics  : channel.empty()
        pathseq_score_warnings = run_pathseq ? PATHSEQ.out.score_warnings : channel.empty()
        mixcr_clns     = run_mixcr   ? MIXCR.out.clns               : channel.empty()
        mixcr_report   = run_mixcr   ? MIXCR.out.report             : channel.empty()
        mixcr_report_json = run_mixcr ? MIXCR.out.report_json       : channel.empty()
        mixcr_step_reports_txt  = run_mixcr ? MIXCR.out.step_reports_txt  : channel.empty()
        mixcr_step_reports_json = run_mixcr ? MIXCR.out.step_reports_json : channel.empty()
        mixcr_clonotypes = run_mixcr ? MIXCR_EXPORT_CLONES.out.clonotypes_tsv : channel.empty()
        telseq         = run_telseq  ? TELSEQ.out.telseq            : channel.empty()
}
