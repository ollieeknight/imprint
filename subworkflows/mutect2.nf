include {
    MUTECT2_CALL
    MUTECT2_PILEUP
    MUTECT2_GATHER_PILEUPS
    MUTECT2_CONTAMINATION
    MUTECT2_MERGE_STATS
    MUTECT2_MERGE_VCFS
    MUTECT2_FILTER
} from '../modules/mutect2'
include { fastaRef; indexed; tumors; normalsOnce } from './common'

// Shards finish in any order; each group closes once all of them have arrived.
def unshard(meta, Map extra = [:]) {
    groupKey([meta: meta.subMap(meta.keySet() - 'shard_count')] + extra, meta.shard_count)
}

workflow MUTECT2 {
    take:
        ch_paired_bams // [meta, tumor_bam, tumor_bai, normal_bam, normal_bai]
        ch_shards      // [interval_gz, interval_tbi], one per scatter shard

    main:
        def fasta = fastaRef()

        def ch_shard_count = ch_shards.count()

        def ch_scattered = ch_paired_bams.combine(ch_shards).combine(ch_shard_count)

        MUTECT2_CALL(
            ch_scattered.map { meta, tb, tbai, nb, nbai, gz, tbi, n -> [meta + [shard_count: n], tb, tbai, nb, nbai, gz, tbi] },
            fasta, indexed(params.gnomad_germline_resource_vcf), indexed(params.pon_vcf)
        )

        def ch_vcfs = MUTECT2_CALL.out.unfiltered_vcf
            .map { meta, vcf, tbi -> [unshard(meta), vcf, tbi] }
            .groupTuple(by: 0)
            .map { key, vcfs, tbis -> [key.target.meta, vcfs.sort { f -> f.name }, tbis.sort { f -> f.name }] }
        MUTECT2_MERGE_VCFS(ch_vcfs)

        def ch_f1r2 = MUTECT2_CALL.out.f1r2
            .map { meta, f1r2 -> [unshard(meta), f1r2] }
            .groupTuple(by: 0)
            .map { key, f1r2s -> [key.target.meta, f1r2s.sort { f -> f.name }] }

        MUTECT2_MERGE_STATS(
            MUTECT2_CALL.out.stats
                .map { meta, stats -> [unshard(meta), stats] }
                .groupTuple(by: 0)
                .map { key, stats -> [key.target.meta, stats.sort { f -> f.name }] }
        )

        // Normals are piled up once per shard, not once per pair.
        MUTECT2_PILEUP(
            tumors(ch_paired_bams).mix(normalsOnce(ch_paired_bams))
                .combine(ch_shards)
                .combine(ch_shard_count)
                .map { meta, role, sample_id, bam, bai, gz, tbi, n -> [meta + [shard_count: n], role, sample_id, bam, bai, gz, tbi] },
            indexed(params.gnomad_pileup_summaries_vcf)
        )

        MUTECT2_GATHER_PILEUPS(
            MUTECT2_PILEUP.out.pileup
                .map { meta, role, sample_id, table -> [unshard(meta, [role: role, sample_id: sample_id]), table] }
                .groupTuple(by: 0)
                .map { key, tables -> [key.target.meta, key.target.role, key.target.sample_id, tables.sort { f -> f.name }] },
            fasta
        )

        def ch_pileups = MUTECT2_GATHER_PILEUPS.out.pileup.branch { _meta, role, _sample_id, _table ->
            tumor: role == 'tumor'
            normal: role == 'normal'
        }
        MUTECT2_CONTAMINATION(
            ch_pileups.tumor
                .map { meta, _role, _sample_id, table -> [meta.normal_id, meta, table] }
                .combine(ch_pileups.normal.map { _meta, _role, normal_id, table -> [normal_id, table] }, by: 0)
                .map { _normal_id, meta, tumor_table, normal_table -> [meta, tumor_table, normal_table] }
        )

        def ch_filter_input = MUTECT2_MERGE_VCFS.out.vcf
            .map { meta, vcf, tbi -> [meta.pair_id, meta, vcf, tbi] }
            .join(ch_f1r2.map { meta, f1r2s -> [meta.pair_id, f1r2s] }, failOnDuplicate: true, failOnMismatch: true)
            .join(MUTECT2_CONTAMINATION.out.contamination.map { meta, table -> [meta.pair_id, table] }, failOnDuplicate: true, failOnMismatch: true)
            .join(MUTECT2_CONTAMINATION.out.segments.map { meta, table -> [meta.pair_id, table] }, failOnDuplicate: true, failOnMismatch: true)
            .join(MUTECT2_MERGE_STATS.out.stats.map { meta, stats -> [meta.pair_id, stats] }, failOnDuplicate: true, failOnMismatch: true)
            .map { _pair_id, meta, vcf, tbi, f1r2s, contamination, segments, stats -> [meta, vcf, tbi, f1r2s, contamination, segments, stats] }
        MUTECT2_FILTER(ch_filter_input, fasta)

    emit:
        vcf           = MUTECT2_FILTER.out.vcf
        raw_vcf       = MUTECT2_MERGE_VCFS.out.vcf
        contamination = MUTECT2_CONTAMINATION.out.contamination
        segments      = MUTECT2_CONTAMINATION.out.segments
}
