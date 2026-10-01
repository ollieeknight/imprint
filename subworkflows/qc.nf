include { MOSDEPTH; VERIFYBAMID2; SOMALIER_EXTRACT; SOMALIER_RELATE; RIKER_QC; PER_BASE_ERROR_RATE; MULTIQC } from '../modules/qc'
include { fastaRef; indexed; optionalFile; callingBed } from './common'

workflow QC {
    take:
        ch_bams          // [meta, bam, bai], analysis BAMs
        ch_align_reports // fastp JSON/HTML and duplicate metrics, from ALIGN.out.reports
        ch_config_yaml   // pipeline config YAML for MultiQC

    main:
        def byName = { a, b -> a.name <=> b.name }
        def fasta = fastaRef()

        MOSDEPTH(ch_bams, callingBed())
        def ch_error_metrics = channel.empty()
        if (!params.dupcaller) {
            PER_BASE_ERROR_RATE(ch_bams, fasta, indexed(params.dbsnp))
            ch_error_metrics = PER_BASE_ERROR_RATE.out.metrics
        }
        VERIFYBAMID2(ch_bams, fasta)
        RIKER_QC(ch_bams, fasta, [optionalFile(params.bait_intervals), optionalFile(params.target_intervals)])
        SOMALIER_EXTRACT(ch_bams, fasta, indexed(params.somalier_sites))

        // Sorted so the groups file is identical across runs.
        def ch_groups = ch_bams
            .map { meta, _bam, _bai -> [meta.donor, meta.id] }
            .groupTuple()
            .map { _donor, ids -> ids.sort().join(',') }
            .collectFile(name: 'groups.txt', newLine: true, sort: true)

        SOMALIER_RELATE(SOMALIER_EXTRACT.out.extracted.collect(sort: byName), ch_groups)

        def ch_multiqc_files = ch_align_reports
            .mix(MOSDEPTH.out.cov.map { _meta, f -> f })
            .mix(RIKER_QC.out.metrics.map { _meta, f -> f })
            .mix(VERIFYBAMID2.out.selfsm.map { _meta, f -> f })
            .mix(SOMALIER_RELATE.out.results)
            .mix(ch_config_yaml)
            .flatten()
            .collect(sort: byName)

        MULTIQC(ch_multiqc_files, file(params.multiqc_config))

    emit:
        mosdepth_cov    = MOSDEPTH.out.cov
        selfsm          = VERIFYBAMID2.out.selfsm
        riker_metrics   = RIKER_QC.out.metrics
        riker_charts    = RIKER_QC.out.charts
        error_metrics   = ch_error_metrics
        somalier        = SOMALIER_RELATE.out.results
        multiqc_report  = MULTIQC.out.report
        multiqc_data    = MULTIQC.out.data
}
