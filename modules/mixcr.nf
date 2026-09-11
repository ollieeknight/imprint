process MIXCR {
    label 'process_high'
    label 'process_high_cpu'
    label 'process_long'
    tag "${meta.id}"
    container "${params.container_mixcr}"
    publishDir {
        "${params.outdir}/${meta.donor}/samples/${meta.cell_type}/mixcr"
    }, mode: 'copy'

    input:
    tuple val(meta), path(r1s), path(r2s)
    path license_file

    output:
    tuple val(meta), path("${meta.id}.qc.txt"),                          emit: report
    tuple val(meta), path("${meta.id}.qc.json"),                         emit: report_json
    tuple val(meta), path("${meta.id}.contigs.clns"),                    emit: clns
    tuple val(meta), path("${meta.id}.*.report.txt"),                    emit: step_reports_txt
    tuple val(meta), path("${meta.id}.*.report.json"), optional: true,   emit: step_reports_json

    script:
    def r1_list = r1s instanceof List ? r1s : [r1s]
    def r2_list = r2s instanceof List ? r2s : [r2s]
    def r1_files = r1_list.join(' ')
    def r2_files = r2_list.join(' ')
    """
    export MI_LICENSE_FILE="\$(realpath "${license_file}")"

    cat ${r1_files} > combined_R1.fastq.gz
    cat ${r2_files} > combined_R2.fastq.gz

    mixcr analyze exome-seq \\
        --species hsa \\
        --assemble-longest-contigs \\
        -t ${task.cpus} \\
        combined_R1.fastq.gz combined_R2.fastq.gz \\
        "${meta.id}"
    """
}

process MIXCR_EXPORT_CLONES {
    label 'process_medium'
    label 'process_short'
    tag "${meta.id}"
    container "${params.container_mixcr}"
    publishDir {
        "${params.outdir}/${meta.donor}/samples/${meta.cell_type}/mixcr"
    }, mode: 'copy'

    input:
    tuple val(meta), path(clns)
    path license_file

    output:
    tuple val(meta), path("${meta.id}.clonotypes.*.tsv"), optional: true, emit: clonotypes_tsv

    script:
    """
    export MI_LICENSE_FILE="\$(realpath "${license_file}")"

    mixcr exportClones \\
        --filter-out-of-frames \\
        --filter-stops \\
        --split-files-by chain \\
        --not-covered-as-empty \\
        -cloneId \\
        -readCount \\
        -readFraction \\
        -chains \\
        -vHit -dHit -jHit -cHit \\
        -vGene -dGene -jGene -cGene \\
        -vFamily -jFamily \\
        -nFeature CDR3 \\
        -aaFeature CDR3 \\
        -nLength CDR3 \\
        -isProductive CDR3 \\
        -vBestIdentityPercent -jBestIdentityPercent \\
        -nMutationsCount VRegion \\
        "${clns}" \\
        "${meta.id}.clonotypes.tsv"
    """
}
