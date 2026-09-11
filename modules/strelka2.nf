process MANTA {
    label 'process_high'
    label 'process_high_cpu'
    label 'process_very_long'
    tag "${meta.pair_id}"
    container "${params.container_manta}"
    publishDir { "${params.outdir}/${meta.donor}/pairs/${meta.pair_dir}/variant_calling" }, mode: 'copy',
    saveAs: { fn -> fn.startsWith("${meta.pair_id}.manta.") ? fn : null }

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    tuple path(intervals_gz), path(intervals_tbi)

    output:
    tuple val(meta), path("manta/results/variants/candidateSmallIndels.vcf.gz"), path("manta/results/variants/candidateSmallIndels.vcf.gz.tbi"), emit: indels
    tuple val(meta), path("${meta.pair_id}.manta.somaticSV.vcf.gz"), path("${meta.pair_id}.manta.somaticSV.vcf.gz.tbi"), emit: svs

    script:
    def exome_flag       = params.genome ? '' : '--exome'
    def call_regions_flag = params.genome ? '' : "--callRegions \"${intervals_gz}\""
    """
    configManta.py \\
        --normalBam "${normal_bam}" \\
        --tumorBam  "${tumor_bam}" \\
        --referenceFasta "${params.ref_fasta}" \\
        ${exome_flag} \\
        ${call_regions_flag} \\
        --runDir manta

    ./manta/runWorkflow.py -m local -j ${task.cpus}

    mv manta/results/variants/somaticSV.vcf.gz "${meta.pair_id}.manta.somaticSV.vcf.gz"
    mv manta/results/variants/somaticSV.vcf.gz.tbi "${meta.pair_id}.manta.somaticSV.vcf.gz.tbi"
    rm -rf manta/workspace
    """
}

process CALL {
    label 'process_high'
    label 'process_high_cpu'
    label 'process_very_long'
    tag "${meta.pair_id}"
    container "${params.container_strelka}"
    publishDir { "${params.outdir}/${meta.donor}/pairs/${meta.pair_dir}/variant_calling/raw" }, mode: 'copy'

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), path(manta_indels), path(manta_indels_tbi)
    tuple path(intervals_gz), path(intervals_tbi)

    output:
    tuple val(meta), path("${meta.pair_id}.strelka2.unfiltered.snvs.vcf.gz"),   path("${meta.pair_id}.strelka2.unfiltered.snvs.vcf.gz.tbi"),   emit: vcf
    tuple val(meta), path("${meta.pair_id}.strelka2.unfiltered.indels.vcf.gz"), path("${meta.pair_id}.strelka2.unfiltered.indels.vcf.gz.tbi"), emit: indels

    script:
    def exome_flag        = params.genome ? '' : '--exome'
    def call_regions_flag = params.genome ? '' : "--callRegions \"${intervals_gz}\""
    """
    configureStrelkaSomaticWorkflow.py \\
        --normalBam "${normal_bam}" \\
        --tumorBam  "${tumor_bam}" \\
        --referenceFasta "${params.ref_fasta}" \\
        ${call_regions_flag} \\
        --indelCandidates "${manta_indels}" \\
        ${exome_flag} \\
        --runDir strelka2

    sed -i 's/isEmail = isLocalSmtp()/isEmail = False/g' strelka2/runWorkflow.py

    ./strelka2/runWorkflow.py -m local -j ${task.cpus}

    mv strelka2/results/variants/somatic.snvs.vcf.gz "${meta.pair_id}.strelka2.unfiltered.snvs.vcf.gz"
    mv strelka2/results/variants/somatic.snvs.vcf.gz.tbi "${meta.pair_id}.strelka2.unfiltered.snvs.vcf.gz.tbi"
    mv strelka2/results/variants/somatic.indels.vcf.gz "${meta.pair_id}.strelka2.unfiltered.indels.vcf.gz"
    mv strelka2/results/variants/somatic.indels.vcf.gz.tbi "${meta.pair_id}.strelka2.unfiltered.indels.vcf.gz.tbi"
    rm -rf strelka2/workspace
    """
}
