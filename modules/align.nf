process BWA_MEM3_LANE_BULK {
    tag "${meta.id}"
    container "${params.container_align}"
    memory 48.GB
    cpus   { (reads1.size() + reads2.size()) < 5.GB ? 8 : (reads1.size() + reads2.size()) < 15.GB ? 12 : (reads1.size() + reads2.size()) < 60.GB ? 16 : 24 }

    input:
    tuple val(meta), path(reads1), path(reads2)
    tuple val(bwa_prefix), path(bwa_index)

    output:
    tuple val(meta), path("${reads1.name.replace('_R1_trimmed.fastq.gz', '.fixmate.bam')}"), emit: bam

    script:
    def fixmate_bam = reads1.name.replace('_R1_trimmed.fastq.gz', '.fixmate.bam')
    def lane_id     = reads1.name.replace('_R1_trimmed.fastq.gz', '')
    """
    bwa-mem3 mem \\
        -t ${task.cpus} --bam=0 \\
        -Y -K 100000000 \\
        -R "@RG\\tID:${lane_id}\\tSM:${meta.id}\\tPL:ILLUMINA\\tLB:${meta.id}" \\
        "${bwa_prefix}" \\
        ${reads1} ${reads2} \\
    | samtools fixmate -@ ${task.cpus} -m -O bam,level=1 - "${fixmate_bam}"
    """
}

process SORT_LANE_BULK {
    tag "${meta.id}"
    container "${params.container_samtools}"
    memory { bam.size() < 20.GB ? 24.GB : 48.GB }
    cpus   8

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${bam.name.replace('.fixmate.bam', '.bam')}"), emit: bam

    script:
    def lane_bam    = bam.name.replace('.fixmate.bam', '.bam')
    def total_gb    = task.memory ? task.memory.toGiga() : 24
    def sort_mem_gb = Math.max(1, ((total_gb - 4) / task.cpus).intValue())
    """
    samtools sort -@ ${task.cpus} -m ${sort_mem_gb}G -O bam -o "${lane_bam}" ${bam}
    """
}

process MERGE_TAGGED_BAMS {
    label 'process_medium'
    tag "${meta.id}"
    container "${params.container_samtools}"

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("${meta.id}_merged_tagged.bam"), path("${meta.id}_merged_tagged.bam.bai"), emit: tagged_bam

    script:
    def bam_list = bams instanceof List ? bams : [bams]
    def merge_cmd = bam_list.size() == 1
        ? "ln -s \"\$(readlink -f ${bam_list[0]})\" ${meta.id}_merged_tagged.bam"
        : "samtools merge -f -@ ${task.cpus} ${meta.id}_merged_tagged.bam ${bam_list.join(' ')}"
    """
    ${merge_cmd}
    samtools index ${meta.id}_merged_tagged.bam
    """
}

process MARK_DUPLICATES {
    label 'process_high'
    tag "${meta.id}"
    container "${params.container_gatk}"
    publishDir { "${params.outdir}/${meta.donor}/samples/${meta.cell_type}/qc/markdup" }, mode: 'copy', pattern: '*_markdup_metrics.txt'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}_dedup.bam"), path("${meta.id}_dedup.bam.bai"), emit: dedup_bam
    tuple val(meta), path("${meta.id}_markdup_metrics.txt"), emit: metrics

    script:
    def mem_mb = task.memory ? (task.memory.mega * 0.8).intValue() : 12288
    """
    gatk --java-options "-Xmx${mem_mb}M" MarkDuplicates \\
        --INPUT "${bam}" \\
        --OUTPUT "${meta.id}_dedup.bam" \\
        --METRICS_FILE "${meta.id}_markdup_metrics.txt" \\
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE ${params.optical_dup_dist} \\
        --CREATE_INDEX true \\
        --VALIDATION_STRINGENCY SILENT

    mv "${meta.id}_dedup.bai" "${meta.id}_dedup.bam.bai"
    """
}

process EXPORT_CRAM {
    label 'process_high'
    tag "${meta.id}"
    container "${params.container_samtools}"
    publishDir { "${params.outdir}/${meta.donor}/samples/${meta.cell_type}/alignment" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    tuple path(fasta), path(fai), path(dict)

    output:
    tuple val(meta), path("${meta.id}.cram"), path("${meta.id}.cram.crai"), emit: cram

    script:
    """
    samtools view -@ ${task.cpus} -C --output-fmt-option version=3.0 -T "${fasta}" -o "${meta.id}.cram" "${bam}"
    samtools index "${meta.id}.cram"
    """
}

process CORRECT_OVERLAPPING_BASES {
    label 'process_low'
    tag "${meta.id}"
    container "${params.container_align}"
    publishDir { "${params.outdir}/${meta.donor}/samples/${meta.cell_type}/qc/fgbio" }, mode: 'copy', pattern: '*_overlap_metrics.txt'

    input:
    tuple val(meta), path(bam), path(bai)
    tuple path(fasta), path(fai), path(dict)

    output:
    tuple val(meta), path("${meta.id}_clipped.bam"), path("${meta.id}_clipped.bam.bai"), emit: bam
    tuple val(meta), path("${meta.id}_overlap_metrics.txt"), emit: metrics

    script:
    def total_gb = task.memory ? task.memory.toGiga() : 8
    def fgbio_gb = Math.max(4, (total_gb / 2).toInteger())
    def sort_mem_gb = Math.max(1, total_gb.intdiv(task.cpus * 4))
    """
    samtools sort -@ ${task.cpus} -m ${sort_mem_gb}G -n -O bam -o qname_sorted.bam "${bam}"
    fgbio -Xmx${fgbio_gb}g --compression=1 --async-io=true CallOverlappingConsensusBases \\
        --input qname_sorted.bam --output clipped.qname.bam \\
        --metrics "${meta.id}_overlap_metrics.txt" --ref "${fasta}" \\
        --threads ${task.cpus} --agreement-strategy MaxQual --disagreement-strategy MaskLowerQual
    samtools sort -@ ${task.cpus} -m ${sort_mem_gb}G -O bam -o "${meta.id}_clipped.bam" clipped.qname.bam
    samtools index "${meta.id}_clipped.bam"
    """
}
