process MUTECT2_NORMAL_ONLY {
    label 'process_medium'
    label 'process_high_memory'
    label 'process_long'
    tag "${meta.id}"
    container "${params.container_gatk}"
    publishDir "${params.outdir}/cohort/pon/normals", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.pon.vcf.gz"), path("${meta.id}.pon.vcf.gz.tbi"), emit: vcf

    script:
    def mem_gb = task.memory ? task.memory.toGiga() - 4 : 8
    def effective_bed = params.off_target ? params.padded_intervals_bed : params.intervals_bed
    def intervals_flag = params.genome ? '' : "-L \"${effective_bed}\""
    """
    gatk --java-options "-Xmx${mem_gb}g" Mutect2 \\
        -R "${params.ref_fasta}" \\
        -I "${bam}" \\
        ${intervals_flag} \\
        --germline-resource "${params.gnomad_germline_resource_vcf}" \\
        --max-mnp-distance 0 \\
        --native-pair-hmm-threads ${task.cpus} \\
        -O "${meta.id}.pon.vcf.gz"
    """
}

process GENOMICSDB_IMPORT_PON {
    label 'process_medium'
    label 'process_high_memory'
    tag 'cohort'
    container "${params.container_gatk}"

    input:
    tuple val(cohort), path(vcfs), path(tbis)

    output:
    tuple val(cohort), path("pon_db"), emit: genomicsdb

    script:
    def mem_gb = task.memory ? task.memory.toGiga() - 4 : 16
    def effective_bed = params.off_target ? params.padded_intervals_bed : params.intervals_bed
    def vcf_inputs = vcfs.collect { v -> "-V ${v}" }.join(" \\\n        ")
    // GenomicsDBImport always needs intervals. A sequence dictionary is not a
    // valid -L argument, so genome mode builds a GATK .intervals list from the
    // FASTA index instead.
    def intervals_flag = params.genome ? '-L genome.intervals' : "-L \"${effective_bed}\""
    def genome_intervals = params.genome
        ? "awk '{ print \$1\":1-\"\$2 }' \"${params.genome_fai}\" > genome.intervals"
        : ''
    """
    ${genome_intervals}

    gatk --java-options "-Xmx${mem_gb}g" GenomicsDBImport \\
        -R "${params.ref_fasta}" \\
        ${intervals_flag} \\
        --genomicsdb-workspace-path pon_db \\
        --merge-input-intervals \\
        ${vcf_inputs}
    """
}

process CREATE_PON {
    label 'process_medium'
    label 'process_high_memory'
    tag 'cohort'
    container "${params.container_gatk}"
    publishDir "${params.outdir}/cohort/pon", mode: 'copy'

    input:
    tuple val(cohort), path(pon_db)

    output:
    tuple val(cohort), path("pon.vcf.gz"), path("pon.vcf.gz.tbi"), emit: pon

    script:
    def mem_gb = task.memory ? task.memory.toGiga() - 4 : 16
    """
    gatk --java-options "-Xmx${mem_gb}g" CreateSomaticPanelOfNormals \\
        -R "${params.ref_fasta}" \\
        --min-sample-count 2 \\
        --germline-resource "${params.gnomad_germline_resource_vcf}" \\
        -V "gendb://${pon_db}" \\
        -O pon.vcf.gz

    tabix -f -p vcf pon.vcf.gz
    """
}
