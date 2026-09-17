include {
    CALL
    PILEUP
    GATHER_PILEUPS
    CONTAMINATION
    MERGE_STATS
    FILTER
    MERGE_VCFS
} from '../modules/mutect2'

workflow MUTECT2 {
    take:
        ch_paired_bams_scattered // [meta+{interval_count}, tb, tbai, nb, nbai, interval_gz, interval_tbi]

    main:
        CALL(ch_paired_bams_scattered)

        CALL.out.unfiltered_vcf
            .map { meta, vcf, tbi ->
                def clean_meta = meta.subMap(meta.keySet() - 'interval_count')
                [groupKey(clean_meta, meta.interval_count), vcf, tbi]
            }
            .groupTuple(by: 0)
            .map { gkey, vcfs, tbis ->
                def sorted = [vcfs, tbis].transpose().sort { a, b -> a[0].name <=> b[0].name }.transpose()
                [gkey.target, 'mutect2', sorted[0], sorted[1]]
            }
            .set { ch_mutect2_vcfs_to_merge }

        MERGE_VCFS(ch_mutect2_vcfs_to_merge)

        CALL.out.f1r2
            .map { meta, f1r2 ->
                def clean_meta = meta.subMap(meta.keySet() - 'interval_count')
                [groupKey(clean_meta, meta.interval_count), f1r2]
            }
            .groupTuple(by: 0)
            .map { gkey, f1r2s -> [gkey.target, f1r2s.sort { f -> f.name }] }
            .set { ch_mutect2_grouped_f1r2 }

        CALL.out.stats
            .map { meta, stats ->
                def clean_meta = meta.subMap(meta.keySet() - 'interval_count')
                [groupKey(clean_meta, meta.interval_count), stats]
            }
            .groupTuple(by: 0)
            .map { gkey, stats -> [gkey.target, stats.flatten().sort { f -> f.name }] }
            .set { ch_mutect2_grouped_stats }

        MERGE_STATS(ch_mutect2_grouped_stats)

        ch_paired_bams_scattered.multiMap { meta, tb, tbai, nb, nbai, interval_gz, interval_tbi ->
            tumour: [meta, 'tumour', meta.tumor_id,  tb, tbai, interval_gz, interval_tbi]
            normal: [meta, 'normal', meta.normal_id, nb, nbai, interval_gz, interval_tbi]
        }.set { ch_pileup_scattered }

        PILEUP(ch_pileup_scattered.tumour.mix(ch_pileup_scattered.normal))

        PILEUP.out.pileup
            .map { meta, role, sample_id, shard ->
                def clean_meta = meta.subMap(meta.keySet() - 'interval_count')
                [groupKey([meta: clean_meta, role: role, sample_id: sample_id], meta.interval_count), shard]
            }
            .groupTuple(by: 0)
            .map { gkey, shards ->
                [gkey.target.meta, gkey.target.role, gkey.target.sample_id, shards.sort { f -> f.name }]
            }
            .set { ch_pileups_to_gather }

        GATHER_PILEUPS(ch_pileups_to_gather)

        GATHER_PILEUPS.out.pileup
            .branch { _meta, role, _sample_id, _pileup ->
                tumour: role == 'tumour'
                normal: role == 'normal'
            }
            .set { ch_pileups }

        ch_pileups.tumour
            .map { meta, _role, _sample_id, pileup -> [meta.pair_id, meta, pileup] }
            .join(
                ch_pileups.normal.map { meta, _role, _sample_id, pileup -> [meta.pair_id, pileup] },
                by: 0, failOnDuplicate: true, failOnMismatch: true
            )
            .map { _pair_id, meta, t_pileup, n_pileup -> [meta, t_pileup, n_pileup] }
            .set { ch_contamination_input }

        CONTAMINATION(ch_contamination_input)

        def ch_contamination_split = CONTAMINATION.out.contamination.multiMap { meta, table ->
            filtering: [meta, table]
            manifest:  [meta, table]
        }
        def ch_segments_split = CONTAMINATION.out.segments.multiMap { meta, table ->
            filtering: [meta, table]
            manifest:  [meta, table]
        }
        def ch_merged_vcf_split = MERGE_VCFS.out.vcf.multiMap { meta, vcf, tbi ->
            filtering: [meta, vcf, tbi]
            manifest: [meta, vcf, tbi]
        }

        ch_merged_vcf_split.filtering
            .map { meta, vcf, tbi -> [meta.pair_id, meta, vcf, tbi] }
            .join(ch_mutect2_grouped_f1r2.map { meta, f1r2s -> [meta.pair_id, f1r2s] },
                by: 0, failOnDuplicate: true, failOnMismatch: true)
            .join(ch_contamination_split.filtering.map { meta, c -> [meta.pair_id, c] },
                by: 0, failOnDuplicate: true, failOnMismatch: true)
            .join(ch_segments_split.filtering.map { meta, s -> [meta.pair_id, s] },
                by: 0, failOnDuplicate: true, failOnMismatch: true)
            .join(MERGE_STATS.out.stats.map { meta, s -> [meta.pair_id, s] },
                by: 0, failOnDuplicate: true, failOnMismatch: true)
            .map { _pair_id, meta, vcf, tbi, f1r2s, contamination, segments, stats ->
                [meta, vcf, tbi, f1r2s, contamination, segments, stats]
            }
            .set { ch_filter_inputs }

        FILTER(ch_filter_inputs)

    emit:
        vcf           = FILTER.out.vcf
        raw_vcf       = ch_merged_vcf_split.manifest
        contamination = ch_contamination_split.manifest
        segments      = ch_segments_split.manifest
}
