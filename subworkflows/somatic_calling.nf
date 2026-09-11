include { PREPARE_INTERVALS; SPLIT_INTERVALS } from '../modules/intervals'
include { MUTECT2 }                               from './mutect2'
include { STRELKA2 }                              from './strelka2'
include { DEEPSOMATIC }                           from './deepsomatic'
include { BASE_RECALIBRATOR; APPLY_BQSR }         from '../modules/mutect2'

workflow SOMATIC_CALLING {
    take:
        ch_paired_bams // [meta, tumor_bam, tumor_bai, normal_bam, normal_bai]

    main:
        // 0. Prepare compressed intervals for Manta and Strelka2
        // In WGS mode: Manta/Strelka2 run genome-wide (no --callRegions); provide a dummy
        // channel so downstream processes still receive two path values (scripts check params.genome).
        // In WES mode: compress the calling BED (padded in off-target mode, standard otherwise).
        def effective_bed = params.off_target ? params.padded_intervals_bed : params.intervals_bed
        if (params.dupcaller) {
            ch_intervals_gz = channel.empty()
            ch_intervals_bed_split = channel.empty()
        } else if (params.genome) {
            ch_intervals_gz = channel.value([file('NO_FILE_GZ'), file('NO_FILE_TBI')])
            ch_intervals_bed_split = channel.value(file('NO_FILE'))
        } else {
            ch_intervals_bed = channel.value(file(effective_bed))
            PREPARE_INTERVALS(ch_intervals_bed)
            ch_intervals_gz = PREPARE_INTERVALS.out.intervals
            ch_intervals_bed_split = channel.value(file(effective_bed))
        }

        // 0b. Split intervals for scatter-gather (Mutect2 only)
        // Genome mode uses NO_FILE and genome_scatter_count.
        SPLIT_INTERVALS(ch_intervals_bed_split)
        ch_split_source = SPLIT_INTERVALS.out

        // Collect interval pairs with their count.
        def ch_intervals_with_count = ch_split_source.intervals_gz
            .flatten()
            .map { f -> [f.simpleName, f] }
            .join(
                ch_split_source.intervals_tbi
                    .flatten()
                    .map { f -> [f.simpleName, f] },
                by: 0, failOnDuplicate: true, failOnMismatch: true
            )
            .map { _name, gz, tbi -> [gz, tbi] }
            .collect()
            .flatMap { paths ->
                def pairs = paths.collate(2)
                def count = pairs.size()
                pairs.collect { pair -> [pair[0], pair[1], count] }
            }

        // Use the same recalibrated BAMs for bulk calling and read-evidence annotation.
        if (!params.skip_bqsr) {
            // Tumour: one entry per pair
            ch_paired_bams
                .map { meta, tb, tbai, _nb, _nbai -> [meta, 'tumour', meta.tumor_id, tb, tbai] }
                .set { ch_tumour_for_bqsr }

            // Normal: deduplicate to 1 per donor (BQSR output is pair-independent)
            ch_paired_bams
                .map { meta, _tb, _tbai, nb, nbai -> [meta.normal_id, meta, nb, nbai] }
                .groupTuple(by: 0)
                .map { normal_id, metas, nbams, nbais -> [metas[0], 'normal', normal_id, nbams[0], nbais[0]] }
                .set { ch_normal_for_bqsr }

            // Each tumour pair needs its donor's recalibrated normal BAM.
            ch_normal_for_bqsr.mix(ch_tumour_for_bqsr).set { ch_bams_for_bqsr }

            BASE_RECALIBRATOR(ch_bams_for_bqsr)

            def ch_recal_table_split = BASE_RECALIBRATOR.out.recal_table.multiMap { meta, role, sample_id, table ->
                apply:    [meta, role, sample_id, table]
                manifest: [meta, role, sample_id, table]
            }

            ch_bams_for_bqsr
                .map { meta, role, sample_id, bam, bai -> [ sample_id, meta, role, bam, bai ] }
                .join(
                    ch_recal_table_split.apply
                        .map { _meta, _role, sample_id, tbl -> [ sample_id, tbl ] },
                    by: 0, failOnDuplicate: true, failOnMismatch: true
                )
                .map { sample_id, meta, role, bam, bai, tbl -> [ meta, role, sample_id, bam, bai, tbl ] }
                .set { ch_apply_bqsr_input }

            APPLY_BQSR(ch_apply_bqsr_input)

            APPLY_BQSR.out.bam
                .branch { _meta, role, _sample_id, _bam, _bai ->
                    tumour: role == 'tumour'
                    normal: role == 'normal'
                }
                .set { ch_bqsr_bams }

            // Fan deduplicated normal back to each tumour pair
            ch_bqsr_bams.tumour
                .map { meta, _role, _sample_id, tb, tbai -> [ meta.donor, meta, tb, tbai ] }
                .combine(
                    ch_bqsr_bams.normal
                        .map { meta, _role, _sample_id, nb, nbai -> [ meta.donor, nb, nbai ] },
                    by: 0
                )
                .map { _donor, meta, tb, tbai, nb, nbai -> [ meta, tb, tbai, nb, nbai ] }
                .set { ch_paired_bams_for_mutect2 }

            ch_bqsr_recal = ch_recal_table_split.manifest

        } else {
            ch_paired_bams.set { ch_paired_bams_for_mutect2 }
            ch_bqsr_recal = channel.empty()
        }

        // Fan out recalibrated BAMs, or analysis BAMs when BQSR is off.
        def ch_calling_bams = ch_paired_bams_for_mutect2.multiMap { meta, tb, tbai, nb, nbai ->
            mutect2:     [meta, tb, tbai, nb, nbai]
            strelka2:    [meta, tb, tbai, nb, nbai]
            deepsomatic: [meta, tb, tbai, nb, nbai]
            annotation:  [meta, tb, tbai, nb, nbai]
        }

        ch_paired_bams_mutect2_scattered = ch_calling_bams.mutect2
            .combine(ch_intervals_with_count)
            .map { meta, tb, tbai, nb, nbai, gz, tbi, count ->
                [meta + [interval_count: count], tb, tbai, nb, nbai, gz, tbi]
            }

        MUTECT2(ch_paired_bams_mutect2_scattered)

        STRELKA2(ch_calling_bams.strelka2, ch_intervals_gz)

        // DeepSomatic requires a GPU and contributes to the PASS union.
        DEEPSOMATIC(ch_calling_bams.deepsomatic)

    emit:
        // VAFATOR and VarLociraptor use the same BAMs as the callers.
        calling_bams    = ch_calling_bams.annotation
        mutect2_vcf     = MUTECT2.out.vcf
        mutect2_raw_vcf = MUTECT2.out.raw_vcf
        mutect2_contamination = MUTECT2.out.contamination
        mutect2_segments = MUTECT2.out.segments
        strelka_snv     = STRELKA2.out.snv
        strelka_indel   = STRELKA2.out.indel
        manta_sv        = STRELKA2.out.manta_sv
        deepsomatic_vcf = DEEPSOMATIC.out
        bqsr_recal      = ch_bqsr_recal
}
