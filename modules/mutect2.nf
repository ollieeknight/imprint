process CALL {
    label 'process_low'
    label 'process_high_memory'
    label 'process_long'
    tag "${meta.pair_id}"
    container "${params.container_gatk}"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), path(interval_gz), path(interval_tbi)

    output:
    tuple val(meta), path("${meta.tumor_id}_mutect2.*.vcf.gz"),      path("${meta.tumor_id}_mutect2.*.vcf.gz.tbi"), emit: unfiltered_vcf
    tuple val(meta), path("${meta.tumor_id}_mutect2.*.vcf.gz.stats"),                                              emit: stats
    tuple val(meta), path("${meta.tumor_id}_*_f1r2.tar.gz"),                                                       emit: f1r2

    script:
    def mem_gb            = task.memory ? task.memory.toGiga() - 2 : 14
    def pcr_indel_model   = '--pcr-indel-model AGGRESSIVE'
    // null autodetects short inserts from fastp.
    def auto_short       = (params.filter_soft_clips != false) && (meta.short_inserts ?: false)
    def soft_clip_filter = params.filter_soft_clips == true || auto_short
        ? '--dont-use-soft-clipped-bases' : ''
    def pon_flags = "--panel-of-normals ${params.pon_vcf}"
    """
    gatk --java-options "-Xmx${mem_gb}g" Mutect2 \\
        -R "${params.ref_fasta}" \\
        -I "${tumor_bam}"  --tumor-sample  "${meta.tumor_id}" \\
        -I "${normal_bam}" --normal-sample "${meta.normal_id}" \\
        --germline-resource "${params.gnomad_germline_resource_vcf}" \\
        ${pon_flags} \\
        -L "${interval_gz}" \\
        --f1r2-tar-gz "${meta.tumor_id}_${interval_gz.simpleName}_f1r2.tar.gz" \\
        -O "${meta.tumor_id}_mutect2.${interval_gz.simpleName}.unfiltered.vcf.gz" \\
        ${pcr_indel_model} \\
        ${soft_clip_filter} \\
        --tumor-lod-to-emit 0 \\
        --initial-tumor-lod 0 \\
        --max-mnp-distance 2 \\
        --native-pair-hmm-threads ${task.cpus} \\
        --tmp-dir .
    """
}

process PILEUP {
    label 'process_low'
    tag "${sample_id}"
    container "${params.container_gatk}"

    input:
    tuple val(meta), val(role), val(sample_id), path(bam), path(bai), path(interval_gz), path(interval_tbi)

    output:
    tuple val(meta), val(role), val(sample_id), path("${sample_id}_${interval_gz.simpleName}.pileup.table"), emit: pileup

    script:
    def mem_gb = task.memory ? task.memory.toGiga() - 2 : 6
    """
    gatk --java-options "-Xmx${mem_gb}g" GetPileupSummaries \\
        -I "${bam}" \\
        -V "${params.gnomad_pileup_summaries_vcf}" \\
        -L "${interval_gz}" \\
        -O "${sample_id}_${interval_gz.simpleName}.pileup.table"
    """
}

process GATHER_PILEUPS {
    label 'process_low'
    label 'process_very_short'
    tag "${meta.pair_id}"
    container "${params.container_gatk}"

    input:
    tuple val(meta), val(role), val(sample_id), path(pileup_tables)

    output:
    tuple val(meta), val(role), val(sample_id), path("${sample_id}_pileups.table"), emit: pileup

    script:
    def mem_gb     = task.memory ? task.memory.toGiga() - 2 : 6
    def input_args = pileup_tables.collect { v -> "-I ${v}" }.join(' ')
    """
    gatk --java-options "-Xmx${mem_gb}g" GatherPileupSummaries \\
        ${input_args} \\
        --sequence-dictionary "${params.ref_dict}" \\
        -O "${sample_id}_pileups.table"

    rm ${pileup_tables.join(' ')}
    """
}

process CONTAMINATION {
    label 'process_low'
    label 'process_very_short'
    tag "${meta.pair_id}"
    container "${params.container_gatk}"
    publishDir { "${params.outdir}/${meta.donor}/pairs/${meta.pair_dir}/variant_calling/raw" }, mode: 'copy'

    input:
    tuple val(meta), path(tumor_pileup), path(normal_pileup)

    output:
    tuple val(meta), path("${meta.pair_id}.contamination.table"), emit: contamination
    tuple val(meta), path("${meta.pair_id}.segments.table"),      emit: segments

    script:
    def mem_gb = task.memory ? task.memory.toGiga() - 2 : 6
    """
    gatk --java-options "-Xmx${mem_gb}g" CalculateContamination \\
        -I "${tumor_pileup}" \\
        -matched "${normal_pileup}" \\
        -O "${meta.pair_id}.contamination.table" \\
        --tumor-segmentation "${meta.pair_id}.segments.table"
    """
}

process MERGE_STATS {
    label 'process_low'
    label 'process_low_memory'
    label 'process_very_short'
    tag "${meta.pair_id}"
    container "${params.container_gatk}"

    input:
    tuple val(meta), path(stats)

    output:
    tuple val(meta), path("${meta.tumor_id}_merged.stats"), emit: stats

    script:
    def mem_gb     = task.memory ? task.memory.toGiga() - 2 : 2
    def stats_args = stats.collect { v -> "-stats ${v}" }.join(' ')
    """
    gatk --java-options "-Xmx${mem_gb}g" MergeMutectStats \\
        ${stats_args} \\
        -O "${meta.tumor_id}_merged.stats"

    rm ${stats.join(' ')}
    """
}

process FILTER {
    label 'process_low'
    tag "${meta.pair_id}"
    container "${params.container_gatk}"
    publishDir { "${params.outdir}/${meta.donor}/pairs/${meta.pair_dir}/variant_calling/raw" }, mode: 'copy'

    input:
    tuple val(meta), path(unfiltered_vcf), path(unfiltered_tbi), path(f1r2s), path(contamination), path(segments), path(stats)

    output:
    tuple val(meta), path("${meta.pair_id}.mutect2.filtered.vcf.gz"), path("${meta.pair_id}.mutect2.filtered.vcf.gz.tbi"), emit: vcf

    script:
    def mem_gb        = task.memory ? task.memory.toGiga() - 2 : 6
    def f1r2_inputs   = f1r2s.collect { v -> "-I ${v}" }.join(' ')
    // Keep FilterMutectCalls' native evidence and artefact filters.
    """
    gatk --java-options "-Xmx${mem_gb}g" LearnReadOrientationModel \\
        ${f1r2_inputs} \\
        -O "${meta.tumor_id}_read_orientation_model.tar.gz"

    gatk --java-options "-Xmx${mem_gb}g" FilterMutectCalls \\
        -R "${params.ref_fasta}" \\
        -V "${unfiltered_vcf}" \\
        --stats "${stats}" \\
        --ob-priors "${meta.tumor_id}_read_orientation_model.tar.gz" \\
        --contamination-table "${contamination}" \\
        --tumor-segmentation "${segments}" \\
        -O "${meta.pair_id}.mutect2.filtered.vcf.gz"

    gatk IndexFeatureFile -I "${meta.pair_id}.mutect2.filtered.vcf.gz"
    """
}

process BASE_RECALIBRATOR {
    label 'process_low'
    label 'process_high_memory'
    label 'process_long'
    tag "${sample_id}"
    container "${params.container_gatk}"
    publishDir {
        def cell_type = role == 'tumour' ? meta.tumor_cell_type : meta.normal_cell_type
        "${params.outdir}/${meta.donor}/samples/${cell_type}/qc/bqsr"
    }, mode: 'copy'

    input:
    tuple val(meta), val(role), val(sample_id), path(bam), path(bai)

    output:
    tuple val(meta), val(role), val(sample_id), path("${sample_id}_recal.table"), emit: recal_table

    script:
    def mem_gb        = task.memory ? task.memory.toGiga() - 2 : 8
    def intervals_flag = params.genome ? '' : "-L \"${params.off_target ? params.padded_intervals_bed : params.intervals_bed}\""
    """
    gatk --java-options "-Xmx${mem_gb}g" BaseRecalibrator \\
        -I "${bam}" \\
        -R "${params.ref_fasta}" \\
        ${intervals_flag} \\
        --known-sites "${params.dbsnp}" \\
        --known-sites "${params.known_indels_mills}" \\
        --known-sites "${params.known_snps_1000g}" \\
        -O "${sample_id}_recal.table"

    """
}

process APPLY_BQSR {
    label 'process_low'
    label 'process_high_memory'
    label 'process_long'
    tag "${sample_id}"
    container "${params.container_gatk}"

    input:
    tuple val(meta), val(role), val(sample_id), path(bam), path(bai), path(recal_table)

    output:
    tuple val(meta), val(role), val(sample_id), path("${sample_id}_bqsr.bam"), path("${sample_id}_bqsr.bam.bai"), emit: bam

    script:
    def mem_gb = task.memory ? task.memory.toGiga() - 2 : 8
    """
    gatk --java-options "-Xmx${mem_gb}g" ApplyBQSR \\
        -R "${params.ref_fasta}" \\
        -I "${bam}" \\
        --bqsr-recal-file "${recal_table}" \\
        -O "${sample_id}_bqsr.bam"

    # GATK writes a companion "<prefix>.bai"; the rest of the pipeline addresses indexes
    # as "<file>.bam.bai". Re-index and drop the stray so the emitted pair is consistent.
    rm -f "${sample_id}_bqsr.bai"
    samtools index "${sample_id}_bqsr.bam"
    """
}

process MERGE_VCFS {
    label 'process_low'
    label 'process_low_memory'
    label 'process_very_short'
    tag "${meta.pair_id}"
    container "${params.container_bcftools}"
    publishDir { "${params.outdir}/${meta.donor}/pairs/${meta.pair_dir}/variant_calling/raw" }, mode: 'copy'

    input:
    tuple val(meta), val(caller), path(vcfs), path(tbis)

    output:
    tuple val(meta), path("${meta.pair_id}.${caller}.unfiltered.vcf.gz"), path("${meta.pair_id}.${caller}.unfiltered.vcf.gz.tbi"), emit: vcf

    script:
    def vcf_args = (vcfs instanceof List ? vcfs.sort { v -> v.name } : [vcfs]).join(' ')
    """
    bcftools concat -a -D ${vcf_args} -O z -o "${meta.pair_id}.${caller}.unfiltered.vcf.gz"
    tabix -p vcf "${meta.pair_id}.${caller}.unfiltered.vcf.gz"

    rm ${vcf_args} \$(echo ${vcf_args} | tr ' ' '\n' | sed 's/\$/.tbi/')
    """
}
