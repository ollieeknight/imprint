process FASTP {
    label 'process_medium'
    tag "${meta.id}"
    container "${params.container_fastp}"
    publishDir { "${params.outdir}/${meta.donor}/samples/${meta.cell_type}/qc/fastp" }, mode: 'copy',
        saveAs: { fn -> fn.endsWith('_trimmed.fastq.gz') ? null : fn }

    input:
    tuple val(meta), path(reads1), path(reads2)

    output:
    tuple val(meta), path("${meta.id}_${meta.run_id}_R1_trimmed.fastq.gz"), path("${meta.id}_${meta.run_id}_R2_trimmed.fastq.gz"), optional: true, emit: trimmed_reads
    tuple val(meta), path("${meta.id}_${meta.run_id}_fastp.json"),                                      emit: json
    tuple val(meta), path("${meta.id}_${meta.run_id}_fastp.html"),                                      emit: html

    script:
    def output_args = params.dupcaller ? '' : """--out1 ${meta.id}_${meta.run_id}_R1_trimmed.fastq.gz \\
        --out2 ${meta.id}_${meta.run_id}_R2_trimmed.fastq.gz \\
       """
    def trim_front = params.trim_front as Integer
    def trim_args = trim_front > 0
        ? "--trim_front1 ${trim_front} --trim_front2 ${trim_front} \\\n       "
        : ''
    """
    fastp \\
        --in1 ${reads1} \\
        --in2 ${reads2} \\
        ${output_args}${trim_args} --detect_adapter_for_pe \\
        --length_required 36 \\
        --json ${meta.id}_${meta.run_id}_fastp.json \\
        --html ${meta.id}_${meta.run_id}_fastp.html \\
        --thread ${task.cpus}

    """
}

process MOSDEPTH {
    label 'process_medium'
    label 'process_low_memory'
    tag "${meta.id}"
    container "${params.container_mosdepth}"
    publishDir { "${params.outdir}/${meta.donor}/samples/${meta.cell_type}/qc/mosdepth" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.*"), emit: cov

    script:
    def effective_bed = params.off_target ? params.padded_intervals_bed : params.intervals_bed
    def bed_flag      = params.genome ? '' : "-b ${effective_bed}"
    def quantize      = params.genome ? '0:1:5:10:20:30:50:' : '0:1:5:10:20:50:100:500:'
    """
    mosdepth \\
        -t ${task.cpus} \\
        --fast-mode \\
        --no-per-base \\
        ${bed_flag} \\
        --quantize ${quantize} \\
        ${meta.id} \\
        ${bam}

    """
}

process MULTIQC {
    label 'process_low'
    label 'process_high_memory'
    label 'process_short'
    tag 'cohort'
    container "${params.container_multiqc}"
    publishDir "${params.outdir}/cohort/multiqc", mode: 'copy'

    input:
    path(reports, stageAs: 'reports/?/*')

    output:
    path "multiqc_report.html",  emit: report
    path "multiqc_report_data/", emit: data

    script:
    """
    multiqc reports/ --outdir . --filename multiqc_report.html --config ${params.multiqc_config}

    """
}

process VERIFYBAMID2 {
    label 'process_medium'
    tag "${meta.id}"
    container "${params.container_verifybamid2}"
    publishDir { "${params.outdir}/${meta.donor}/samples/${meta.cell_type}/qc/verifybamid2" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.selfSM"), emit: selfsm

    script:
    def svd_prefix = params.genome ? params.verifybamid2_svd_wgs : params.verifybamid2_svd
    """
    verifybamid2 \\
        --SVDPrefix "${svd_prefix}" \\
        --Reference "${params.ref_fasta}" \\
        --BamFile   "${bam}" \\
        --Output    ${meta.id} \\
        --NumThread ${task.cpus} \\
        --DisableSanityCheck

    """
}

process SOMALIER_EXTRACT {
    label 'process_low'
    label 'process_low_memory'
    label 'process_short'
    tag "${meta.id}"
    container "${params.container_somalier}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    path("${meta.id}.somalier"), emit: extracted

    script:
    """
    somalier extract \\
        -d . \\
        --sites "${params.somalier_sites}" \\
        -f "${params.ref_fasta}" \\
        "${bam}"

    """
}

process SOMALIER_RELATE {
    label 'process_low'
    label 'process_low_memory'
    label 'process_short'
    tag 'cohort'
    container "${params.container_somalier}"
    publishDir "${params.outdir}/cohort/somalier", mode: 'copy'

    input:
    path(extracted_files)
    path(groups_file)

    output:
    path("cohort*"), emit: results

    script:
    """
    somalier relate \\
        --output-prefix cohort \\
        --groups "${groups_file}" \\
        ${extracted_files}
    """
}

process RIKER_QC {
    label 'process_medium'
    tag "${meta.id}"
    container "${params.container_align}"
    publishDir { "${params.outdir}/${meta.donor}/samples/${meta.cell_type}/qc/riker" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai), path(bait_intervals, stageAs: 'baits.bed'), path(target_intervals, stageAs: 'targets.bed')

    output:
    tuple val(meta), path("${meta.id}.*.txt"), emit: metrics
    tuple val(meta), path("${meta.id}.*.pdf"), optional: true, emit: charts

    script:
    if (params.genome) {
        """
        riker multi \\
            -i "${bam}" \\
            -r "${params.ref_fasta}" \\
            -o "${meta.id}" \\
            --tools wgs alignment isize gcbias basic \\
            --threads ${task.cpus}
        """
    } else {
        """
        riker multi \\
            -i "${bam}" \\
            -r "${params.ref_fasta}" \\
            -o "${meta.id}" \\
            --tools hybcap alignment isize gcbias basic \\
            --hybcap::baits "${bait_intervals}" \\
            --hybcap::targets "${target_intervals}" \\
            --threads ${task.cpus}
        """
    }
}

process PER_BASE_ERROR_RATE {
    label 'process_low'
    label 'process_high_memory'
    tag "${meta.id}"
    container "${params.container_align}"
    publishDir { "${params.outdir}/${meta.donor}/samples/${meta.cell_type}/qc/fgbio" },
        mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.error_rate_by_read_position.txt"), emit: metrics

    script:
    def total_gb     = task.memory ? task.memory.toGiga() : 16
    def mem_gb       = Math.max(4, total_gb - 4)
    def variants_arg = params.dbsnp ? "--variants \"${params.dbsnp}\"" : ""
    """
    fgbio -Xmx${mem_gb}g --compression=1 --async-io=true ErrorRateByReadPosition \\
        --input   "${bam}" \\
        --output  "${meta.id}" \\
        --ref     "${params.ref_fasta}" \\
        --min-mapping-quality 20 \\
        ${variants_arg}
    """
}
