process BASE_RECALIBRATOR {
    label 'process_medium'
    tag "${sample_id}"
    container "${params.container_gatk}"
    publishDir {
        def cell_type = role == 'tumor' ? meta.tumor_cell_type : meta.normal_cell_type
        "${params.outdir}/${meta.donor}/samples/${cell_type}/qc/bqsr"
    }, mode: 'copy'

    input:
    tuple val(meta), val(role), val(sample_id), path(bam), path(bai)
    tuple path(fasta), path(fai), path(dict)
    path(calling_bed)
    tuple path(dbsnp), path(dbsnp_tbi), path(mills), path(mills_tbi), path(snps_1000g), path(snps_1000g_tbi)

    output:
    tuple val(meta), val(role), val(sample_id), path("${sample_id}_recal.table"), emit: recal_table

    script:
    def mem_gb         = task.memory ? task.memory.toGiga() - 2 : 8
    def intervals_flag = calling_bed ? "-L \"${calling_bed}\"" : ''
    """
    gatk --java-options "-Xmx${mem_gb}g" BaseRecalibrator \\
        -I "${bam}" \\
        -R "${fasta}" \\
        ${intervals_flag} \\
        --known-sites "${dbsnp}" \\
        --known-sites "${mills}" \\
        --known-sites "${snps_1000g}" \\
        -O "${sample_id}_recal.table"
    """
}

process APPLY_BQSR {
    label 'process_medium'
    tag "${sample_id}"
    container "${params.container_gatk}"

    input:
    tuple val(meta), val(role), val(sample_id), path(bam), path(bai), path(recal_table)
    tuple path(fasta), path(fai), path(dict)

    output:
    tuple val(meta), val(role), val(sample_id), path("${sample_id}_bqsr.bam"), path("${sample_id}_bqsr.bam.bai"), emit: bam

    script:
    def mem_gb = task.memory ? task.memory.toGiga() - 2 : 8
    """
    gatk --java-options "-Xmx${mem_gb}g" ApplyBQSR \\
        -R "${fasta}" \\
        -I "${bam}" \\
        --bqsr-recal-file "${recal_table}" \\
        -O "${sample_id}_bqsr.bam"

    mv "${sample_id}_bqsr.bai" "${sample_id}_bqsr.bam.bai"
    """
}
