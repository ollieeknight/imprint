include { MOSDEPTH; VERIFYBAMID2; SOMALIER_EXTRACT; SOMALIER_RELATE; RIKER_QC; PER_BASE_ERROR_RATE; MULTIQC } from '../modules/qc'
// Riker receives capture intervals only in WES mode.

workflow QC {
    take:
        ch_bams          // [meta, bam, bai], analysis BAMs
        ch_align_reports // fastp JSON/HTML, from ALIGN.out.reports
        ch_config_yaml   // pipeline config YAML for MultiQC

    main:
        def ch_bams_split = ch_bams.multiMap { meta, bam, bai ->
            mosdepth: [meta, bam, bai]
            error:    [meta, bam, bai]
            verify:   [meta, bam, bai]
            riker:    [meta, bam, bai]
            somalier: [meta, bam, bai]
            groups:   [meta, bam, bai]
        }

        MOSDEPTH(ch_bams_split.mosdepth)
        if (params.dupcaller) {
            ch_error_metrics = channel.empty()
        } else {
            PER_BASE_ERROR_RATE(ch_bams_split.error)
            ch_error_metrics = PER_BASE_ERROR_RATE.out.metrics
        }
        VERIFYBAMID2(ch_bams_split.verify)

        RIKER_QC(
            ch_bams_split.riker.map { meta, bam, bai ->
                def baits   = (!params.genome && params.bait_intervals)   ? file(params.bait_intervals)   : file('NO_FILE')
                def targets = (!params.genome && params.target_intervals) ? file(params.target_intervals) : file('NO_FILE')
                [ meta, bam, bai, baits, targets ]
            }
        )

        // Somalier
        SOMALIER_EXTRACT(ch_bams_split.somalier)

        // One donor per line for Somalier relatedness checks.
        def ch_groups = ch_bams_split.groups
            .map { meta, _bam, _bai -> [meta.donor, meta.id] }
            .groupTuple()
            .map { _donor, ids -> ids.join(',') }
            .collectFile(name: 'groups.txt', newLine: true)

        SOMALIER_RELATE(SOMALIER_EXTRACT.out.extracted.collect(), ch_groups)

        def ch_mosdepth_split = MOSDEPTH.out.cov.multiMap { meta, files ->
            multiqc:  [meta, files]
            manifest: [meta, files]
        }

        // Riker metrics feed both Picard and fgbio MultiQC modules.
        def ch_riker_split = RIKER_QC.out.metrics.multiMap { meta, files ->
            multiqc:  [meta, files]
            manifest: [meta, files]
        }

        // Stage configured MultiQC reports.
        ch_align_reports
            .mix(ch_mosdepth_split.multiqc.map { _meta, f -> f })
            .mix(ch_riker_split.multiqc.map { _meta, f -> f })
            .mix(VERIFYBAMID2.out.selfsm.map { _meta, f -> f })
            .mix(SOMALIER_RELATE.out.results)
            .mix(ch_config_yaml)
            .collect()
            .set { ch_multiqc_files }

        MULTIQC(ch_multiqc_files)

    emit:
        mosdepth_cov    = ch_mosdepth_split.manifest
        selfsm          = VERIFYBAMID2.out.selfsm
        riker_metrics   = ch_riker_split.manifest
        riker_charts    = RIKER_QC.out.charts
        error_metrics   = ch_error_metrics
        somalier        = SOMALIER_RELATE.out.results
        multiqc_report  = MULTIQC.out.report
        multiqc_data    = MULTIQC.out.data
}
