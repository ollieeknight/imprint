process PREPARE_INTERVALS {
    label 'process_low'
    tag 'intervals'
    container "${params.container_bcftools}"

    input:
    path(intervals_bed)

    output:
    tuple path("intervals.bed.gz"), path("intervals.bed.gz.tbi"), emit: intervals

    script:
    """
    grep -vE '^(#|track|browser)' "${intervals_bed}" | bgzip -c > intervals.bed.gz
    tabix -p bed intervals.bed.gz
    """
}

process SPLIT_INTERVALS {
    label 'process_low'
    tag 'intervals'
    container "${params.container_gatk}"

    input:
    path(intervals_bed)
    tuple path(fasta), path(fai), path(dict)

    output:
    path("scattered_intervals/*-scattered.bed.gz"),     emit: intervals_gz
    path("scattered_intervals/*-scattered.bed.gz.tbi"), emit: intervals_tbi

    script:
    def scatter_count  = params.genome ? params.genome_scatter_count : params.wes_scatter_count
    def intervals_flag = intervals_bed ? "-L \"${intervals_bed}\"" : ''
    """
    mkdir -p scattered_intervals

    gatk SplitIntervals \\
        -R "${fasta}" \\
        ${intervals_flag} \\
        --scatter-count ${scatter_count} \\
        --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION \\
        -O scattered_intervals

    for f in scattered_intervals/*.interval_list; do
        base=\$(basename \$f .interval_list)
        grep -v '^@' "\$f" \\
            | awk 'BEGIN{OFS="\\t"}{print \$1, \$2-1, \$3}' \\
            > "scattered_intervals/\${base}.bed"
        bgzip "scattered_intervals/\${base}.bed"
        tabix -p bed "scattered_intervals/\${base}.bed.gz"
        rm "\$f"
    done
    """
}
