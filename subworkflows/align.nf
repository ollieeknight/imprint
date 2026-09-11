include { FASTP } from '../modules/qc'

// Remove lane fields before downstream task hashes are calculated.
def sampleMeta(meta) {
    meta.subMap(meta.keySet() - ['run_count', 'run_id'])
}

include {
    BWA_MEM3_LANE_BULK
    SORT_LANE_BULK
    MERGE_TAGGED_BAMS
    MARK_DUPLICATES
    EXPORT_CRAM
    CORRECT_OVERLAPPING_BASES
} from '../modules/align'

include { DUPCALLER_TRIM_LANE; DUPCALLER_ALIGN_LANE; DUPCALLER_SORT_LANE; DUPCALLER_MARK_DUPLICATES; DUPCALLER_VALIDATE_BAM; DUPCALLER_VALIDATE_CRAM } from '../modules/dupcaller'

workflow ALIGN {
    take:
        ch_fastq // [meta, r1, r2], one entry per lane; meta contains run_count

    main:
        def trim_front = params.trim_front as Integer

        // Adapter trimming: adapters detected per pair, reads under 36 bp dropped
        FASTP(ch_fastq)

        // Split fastp outputs for each consumer.
        def ch_fastp_json_split = FASTP.out.json.multiMap { meta, json ->
            short_insert: [meta, json]
            read_length:  [meta, json]
            manifest:     [meta, json]
            report:       [meta, json]
        }
        def ch_fastp_trimmed_split = FASTP.out.trimmed_reads.multiMap { meta, r1, r2 ->
            alignment: [meta, r1, r2]
            mixcr:     [meta, r1, r2]
        }
        def ch_fastp_html_split = FASTP.out.html.multiMap { meta, html ->
            manifest: [meta, html]
            report:   [meta, html]
        }
        if (!params.dupcaller) {
            // Bulk path.
            BWA_MEM3_LANE_BULK(ch_fastp_trimmed_split.alignment)
            SORT_LANE_BULK(BWA_MEM3_LANE_BULK.out.bam)

            // Merge lane BAMs per sample before duplicate marking.
            SORT_LANE_BULK.out.bam
                .map { meta, bam -> [groupKey(meta.id, meta.run_count), meta, bam] }
                .groupTuple(by: 0)
                .map { _id, metas, bams -> [ sampleMeta(metas[0]), bams ] }
                .set { ch_bams_grouped }

            MERGE_TAGGED_BAMS(ch_bams_grouped)

            MARK_DUPLICATES(MERGE_TAGGED_BAMS.out.tagged_bam)

            CORRECT_OVERLAPPING_BASES(MARK_DUPLICATES.out.dedup_bam)
            EXPORT_CRAM(CORRECT_OVERLAPPING_BASES.out.bam)

            ch_analysis_bam  = CORRECT_OVERLAPPING_BASES.out.bam
            ch_merged_bam    = MERGE_TAGGED_BAMS.out.tagged_bam
            ch_markdup_metrics = MARK_DUPLICATES.out.metrics
            ch_overlap_metrics = CORRECT_OVERLAPPING_BASES.out.metrics
            ch_barcode_metrics = channel.empty()
            ch_tag_metrics = channel.empty()
            ch_cram_tag_metrics = channel.empty()
            ch_mode_trimmed_reads = ch_fastp_trimmed_split.mixcr

        } else {
            // xGen UDSeq raw-family path.
            DUPCALLER_TRIM_LANE(
                ch_fastq,
                file(params.dupcaller_umi_allowlist),
                file("${projectDir}/bin/trim_dupcaller_umis.py")
            )
            DUPCALLER_ALIGN_LANE(DUPCALLER_TRIM_LANE.out.reads)
            DUPCALLER_SORT_LANE(DUPCALLER_ALIGN_LANE.out.bam)

            DUPCALLER_SORT_LANE.out.bam
                .map { meta, bam -> [groupKey(meta.id, meta.run_count), meta, bam] }
                .groupTuple(by: 0)
                .map { _id, metas, bams -> [ sampleMeta(metas[0]), bams ] }
                .set { ch_dupcaller_bams_grouped }

            MERGE_TAGGED_BAMS(ch_dupcaller_bams_grouped)
            DUPCALLER_MARK_DUPLICATES(MERGE_TAGGED_BAMS.out.tagged_bam)
            // Tag validation gates downstream DupCaller tasks.
            DUPCALLER_VALIDATE_BAM(
                DUPCALLER_MARK_DUPLICATES.out.bam,
                file(params.dupcaller_umi_allowlist),
                file("${projectDir}/bin/validate_dupcaller_tags.awk")
            )
            EXPORT_CRAM(DUPCALLER_VALIDATE_BAM.out.bam)
            DUPCALLER_VALIDATE_CRAM(
                EXPORT_CRAM.out.cram,
                file(params.dupcaller_umi_allowlist),
                file("${projectDir}/bin/validate_dupcaller_tags.awk")
            )

            ch_analysis_bam  = DUPCALLER_VALIDATE_BAM.out.bam
            ch_merged_bam    = MERGE_TAGGED_BAMS.out.tagged_bam
            ch_markdup_metrics = DUPCALLER_MARK_DUPLICATES.out.metrics
            ch_overlap_metrics = channel.empty()
            ch_barcode_metrics = DUPCALLER_TRIM_LANE.out.metrics
            ch_tag_metrics = DUPCALLER_VALIDATE_BAM.out.metrics
            ch_cram_tag_metrics = DUPCALLER_VALIDATE_CRAM.out.metrics
            ch_mode_trimmed_reads = DUPCALLER_TRIM_LANE.out.reads
        }

    // Detect short inserts from fastp for Mutect2 soft-clip handling.
    ch_short_inserts = ch_fastp_json_split.short_insert
        .map { meta, json -> [groupKey(meta.id, meta.run_count), json] }
        .groupTuple()
        .map { sample_id, jsons ->
            def slurper  = new groovy.json.JsonSlurper()
            def is_short = jsons.any { f ->
                try {
                    def j           = slurper.parse(f)
                    def peak        = (j?.insert_size?.peak ?: 999) as Integer
                    def read_length = ((j?.summary?.before_filtering?.read1_mean_length ?: 150) as Integer) - trim_front
                    peak <= (read_length - 10)
                } catch (_e) { false }
            }
            [sample_id, is_short]
        }

    ch_read_lengths = ch_fastp_json_split.read_length
        .map { meta, json -> [groupKey(meta.id, meta.run_count), json] }
        .groupTuple()
        .map { sample_id, jsons ->
            def slurper = new groovy.json.JsonSlurper()
            def rlen = jsons.collect { f ->
                try {
                    def j = slurper.parse(f)
                    ((j?.summary?.before_filtering?.read1_mean_length ?: 100) as Integer) - trim_front
                } catch (_e) { 100 }
            }.max()
            [sample_id, rlen]
        }

    // fastp_stats: per-sample [id, short_inserts, read_length]
    ch_fastp_stats = ch_short_inserts
        .join(ch_read_lengths)
        .map { id, si, rl -> [id.toString(), (si ?: false), (rl ?: 100)] }

    // Trimmed reads grouped per sample for MiXCR (FASTP runs per lane)
    ch_trimmed_reads_grouped = ch_mode_trimmed_reads
        .map { meta, r1, r2 -> [groupKey(meta.id, meta.run_count), meta, r1, r2] }
        .groupTuple(by: 0)
        .map { _id, metas, r1s, r2s -> [sampleMeta(metas[0]), r1s, r2s] }

    emit:
        analysis_bam   = ch_analysis_bam
        cram           = EXPORT_CRAM.out.cram
        fastp_stats    = ch_fastp_stats
        merged_bam     = ch_merged_bam
        fastp_json     = ch_fastp_json_split.manifest
        fastp_html     = ch_fastp_html_split.manifest
        markdup_metrics = ch_markdup_metrics
        overlap_metrics = ch_overlap_metrics
        barcode_metrics = ch_barcode_metrics
        tag_metrics     = ch_tag_metrics
        cram_tag_metrics = ch_cram_tag_metrics
        trimmed_reads  = ch_trimmed_reads_grouped
        // Include duplicate-marking metrics alongside the fastp reports for MultiQC.
        reports        = ch_fastp_json_split.report.map { _meta, f -> f }
                             .mix(ch_fastp_html_split.report.map { _meta, f -> f })
                             .mix(ch_markdup_metrics.map { _meta, f -> f })
}
