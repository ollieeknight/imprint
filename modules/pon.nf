process MUTECT2_NORMAL_ONLY {
    label 'process_high'
    tag "${meta.id}"
    container "${params.container_gatk}"
    publishDir "${params.outdir}/cohort/pon/normals", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    tuple path(fasta), path(fai), path(dict)
    path(calling_bed)
    tuple path(germline_vcf), path(germline_tbi)

    output:
    tuple val(meta), path("${meta.id}.pon.vcf.gz"), path("${meta.id}.pon.vcf.gz.tbi"), emit: vcf

    script:
    def mem_gb         = task.memory ? task.memory.toGiga() - 4 : 8
    def intervals_flag = calling_bed ? "-L \"${calling_bed}\"" : ''
    """
    gatk --java-options "-Xmx${mem_gb}g" Mutect2 \\
        -R "${fasta}" \\
        -I "${bam}" \\
        ${intervals_flag} \\
        --germline-resource "${germline_vcf}" \\
        --max-mnp-distance 0 \\
        --native-pair-hmm-threads ${task.cpus} \\
        -O "${meta.id}.pon.vcf.gz"
    """
}

process GENOMICSDB_IMPORT_PON {
    label 'process_high'
    tag 'cohort'
    container "${params.container_gatk}"

    input:
    tuple path(vcfs), path(tbis)
    tuple path(fasta), path(fai), path(dict)
    path(calling_bed)

    output:
    path("pon_db"), emit: genomicsdb

    script:
    def mem_gb         = task.memory ? task.memory.toGiga() - 4 : 16
    def vcf_inputs     = vcfs.collect { v -> "-V ${v}" }.join(" \\\n        ")
    def intervals_flag = calling_bed ? "-L \"${calling_bed}\"" : '-L genome.intervals'
    """
    awk '{ print \$1":1-"\$2 }' "${fai}" > genome.intervals

    gatk --java-options "-Xmx${mem_gb}g" GenomicsDBImport \\
        -R "${fasta}" \\
        ${intervals_flag} \\
        --genomicsdb-workspace-path pon_db \\
        --merge-input-intervals \\
        ${vcf_inputs}
    """
}

process CREATE_PON {
    label 'process_high'
    tag 'cohort'
    container "${params.container_gatk}"
    publishDir "${params.outdir}/cohort/pon", mode: 'copy'

    input:
    path(pon_db)
    tuple path(fasta), path(fai), path(dict)
    tuple path(germline_vcf), path(germline_tbi)

    output:
    tuple path("pon.vcf.gz"), path("pon.vcf.gz.tbi"), emit: pon

    script:
    def mem_gb = task.memory ? task.memory.toGiga() - 4 : 16
    """
    gatk --java-options "-Xmx${mem_gb}g" CreateSomaticPanelOfNormals \\
        -R "${fasta}" \\
        --min-sample-count 2 \\
        --germline-resource "${germline_vcf}" \\
        -V "gendb://${pon_db}" \\
        -O pon.vcf.gz

    tabix -f -p vcf pon.vcf.gz
    """
}
