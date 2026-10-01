include { FASTP } from '../modules/qc'

include {
    BWA_MEM3_LANE_BULK
    SORT_LANE_BULK
    MERGE_TAGGED_BAMS
    MARK_DUPLICATES
    EXPORT_CRAM
    CORRECT_OVERLAPPING_BASES
} from '../modules/align'

include { DUPCALLER_TRIM_LANE; DUPCALLER_ALIGN_LANE; DUPCALLER_SORT_LANE; DUPCALLER_MARK_DUPLICATES; DUPCALLER_VALIDATE_BAM; DUPCALLER_VALIDATE_CRAM } from '../modules/dupcaller'

def sampleMeta(meta) {
    meta.subMap(meta.keySet() - ['run_count', 'run_id'])
}

// Lanes finish in any order. Sort them so the merge command, and with it the
// -resume hash of everything downstream, does not change between runs.
def groupLanes(ch) {
    ch.map { meta, bam -> [groupKey(meta.id, meta.run_count), meta, bam] }
        .groupTuple(by: 0)
        .map { _id, metas, bams -> [sampleMeta(metas[0]), bams.sort { b -> b.name }] }
}

def fastpStats(json, int trimFront) {
    try {
        def j    = new groovy.json.JsonSlurper().parse(json)
        def mean = j?.summary?.before_filtering?.read1_mean_length
        def peak = (j?.insert_size?.peak ?: 999) as Integer
        [
            short_insert: peak <= ((mean ?: 150) as Integer) - trimFront - 10,
            read_length:  ((mean ?: 100) as Integer) - trimFront,
        ]
    } catch (_e) {
        [short_insert: false, read_length: 100]
    }
}

workflow ALIGN {
    take:
        ch_fastq

    main:
        def trim_front = params.trim_front as Integer

        FASTP(ch_fastq)

        if (!params.dupcaller) {
            BWA_MEM3_LANE_BULK(FASTP.out.trimmed_reads)
            SORT_LANE_BULK(BWA_MEM3_LANE_BULK.out.bam)
            MERGE_TAGGED_BAMS(groupLanes(SORT_LANE_BULK.out.bam))
            MARK_DUPLICATES(MERGE_TAGGED_BAMS.out.tagged_bam)
            CORRECT_OVERLAPPING_BASES(MARK_DUPLICATES.out.dedup_bam)
            EXPORT_CRAM(CORRECT_OVERLAPPING_BASES.out.bam)

            ch_analysis_bam       = CORRECT_OVERLAPPING_BASES.out.bam
            ch_markdup_metrics    = MARK_DUPLICATES.out.metrics
            ch_overlap_metrics    = CORRECT_OVERLAPPING_BASES.out.metrics
            ch_barcode_metrics    = channel.empty()
            ch_tag_metrics        = channel.empty()
            ch_cram_tag_metrics   = channel.empty()
            ch_mode_trimmed_reads = FASTP.out.trimmed_reads
        } else {
            def allowlist = file(params.dupcaller_umi_allowlist)
            def validator = file("${projectDir}/bin/validate_dupcaller_tags.awk")

            DUPCALLER_TRIM_LANE(ch_fastq, allowlist, file("${projectDir}/bin/trim_dupcaller_umis.py"))
            DUPCALLER_ALIGN_LANE(DUPCALLER_TRIM_LANE.out.reads)
            DUPCALLER_SORT_LANE(DUPCALLER_ALIGN_LANE.out.bam)
            MERGE_TAGGED_BAMS(groupLanes(DUPCALLER_SORT_LANE.out.bam))
            DUPCALLER_MARK_DUPLICATES(MERGE_TAGGED_BAMS.out.tagged_bam)
            // Tag validation gates downstream DupCaller tasks.
            DUPCALLER_VALIDATE_BAM(DUPCALLER_MARK_DUPLICATES.out.bam, allowlist, validator)
            EXPORT_CRAM(DUPCALLER_VALIDATE_BAM.out.bam)
            DUPCALLER_VALIDATE_CRAM(EXPORT_CRAM.out.cram, allowlist, validator)

            ch_analysis_bam       = DUPCALLER_VALIDATE_BAM.out.bam
            ch_markdup_metrics    = DUPCALLER_MARK_DUPLICATES.out.metrics
            ch_overlap_metrics    = channel.empty()
            ch_barcode_metrics    = DUPCALLER_TRIM_LANE.out.metrics
            ch_tag_metrics        = DUPCALLER_VALIDATE_BAM.out.metrics
            ch_cram_tag_metrics   = DUPCALLER_VALIDATE_CRAM.out.metrics
            ch_mode_trimmed_reads = DUPCALLER_TRIM_LANE.out.reads
        }

        // Per sample: [id, short_inserts, read_length]
        ch_fastp_stats = FASTP.out.json
            .map { meta, json -> [groupKey(meta.id, meta.run_count), json] }
            .groupTuple()
            .map { id, jsons ->
                def stats = jsons.collect { json -> fastpStats(json, trim_front) }
                [id.toString(), stats.any { s -> s.short_insert }, stats.collect { s -> s.read_length }.max()]
            }

        // Trimmed reads per sample for MiXCR, in lane order.
        ch_trimmed_reads_grouped = ch_mode_trimmed_reads
            .map { meta, r1, r2 -> [groupKey(meta.id, meta.run_count), meta, r1, r2] }
            .groupTuple(by: 0)
            .map { _id, metas, r1s, r2s ->
                def pairs = [r1s, r2s].transpose().sort { p -> p[0].name }
                [sampleMeta(metas[0]), pairs.collect { p -> p[0] }, pairs.collect { p -> p[1] }]
            }

    emit:
        analysis_bam     = ch_analysis_bam
        cram             = EXPORT_CRAM.out.cram
        fastp_stats      = ch_fastp_stats
        merged_bam       = MERGE_TAGGED_BAMS.out.tagged_bam
        fastp_json       = FASTP.out.json
        fastp_html       = FASTP.out.html
        markdup_metrics  = ch_markdup_metrics
        overlap_metrics  = ch_overlap_metrics
        barcode_metrics  = ch_barcode_metrics
        tag_metrics      = ch_tag_metrics
        cram_tag_metrics = ch_cram_tag_metrics
        trimmed_reads    = ch_trimmed_reads_grouped
        reports          = FASTP.out.json.mix(FASTP.out.html, ch_markdup_metrics).map { _meta, f -> f }
}
