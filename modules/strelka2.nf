process MANTA_PREP_BAM {
    label 'process_medium'
    tag "${sample_id}"
    container "${params.container_samtools}"

    input:
    tuple val(meta), val(role), val(sample_id), path(bam), path(bai)

    output:
    tuple val(meta), val(role), val(sample_id), path("${sample_id}_stripped.bam"), path("${sample_id}_stripped.bam.bai"), emit: bam

    script:
    """
    # Manta 1.6.0 crashes on base qualities above Q70, so cap them at Q70
    # (PHRED+33 'g'). All alignment classes are kept: discordant and
    # supplementary reads are Manta evidence, and each caller filters its own.
    samtools view -h --remove-tag cd,ce,cD,cM,cE \\
        -@ ${task.cpus} \\
        ${bam} \\
        | awk 'BEGIN{OFS="\\t"} /^@/{print;next} {q=\$11; gsub(/[h-~]/,"g",q); \$11=q; print}' \\
        | samtools view -1 -@ ${task.cpus} -o ${sample_id}_stripped.bam
    samtools index -@ ${task.cpus} ${sample_id}_stripped.bam
    """
}

process MANTA {
    label 'process_high_cpu'
    tag "${meta.pair_id}"
    container "${params.container_manta}"
    publishDir { "${params.outdir}/${meta.donor}/pairs/${meta.pair_dir}/variant_calling" }, mode: 'copy',
        saveAs: { fn -> fn.startsWith("${meta.pair_id}.manta.") ? fn : null }

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    tuple path(intervals_gz), path(intervals_tbi)
    tuple path(fasta), path(fai), path(dict)

    output:
    tuple val(meta), path("manta/results/variants/candidateSmallIndels.vcf.gz"), path("manta/results/variants/candidateSmallIndels.vcf.gz.tbi"), emit: indels
    tuple val(meta), path("${meta.pair_id}.manta.somaticSV.vcf.gz"), path("${meta.pair_id}.manta.somaticSV.vcf.gz.tbi"), emit: svs

    script:
    def exome_flag        = params.genome ? '' : '--exome'
    def call_regions_flag = intervals_gz ? "--callRegions \"${intervals_gz}\"" : ''
    """
    configManta.py \\
        --normalBam "${normal_bam}" \\
        --tumorBam  "${tumor_bam}" \\
        --referenceFasta "${fasta}" \\
        ${exome_flag} \\
        ${call_regions_flag} \\
        --runDir manta

    ./manta/runWorkflow.py -m local -j ${task.cpus}

    mv manta/results/variants/somaticSV.vcf.gz "${meta.pair_id}.manta.somaticSV.vcf.gz"
    mv manta/results/variants/somaticSV.vcf.gz.tbi "${meta.pair_id}.manta.somaticSV.vcf.gz.tbi"
    rm -rf manta/workspace
    """
}

process STRELKA2_CALL {
    label 'process_high_cpu'
    tag "${meta.pair_id}"
    container "${params.container_strelka}"
    publishDir { "${params.outdir}/${meta.donor}/pairs/${meta.pair_dir}/variant_calling/raw" }, mode: 'copy'

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), path(manta_indels), path(manta_indels_tbi)
    tuple path(intervals_gz), path(intervals_tbi)
    tuple path(fasta), path(fai), path(dict)

    output:
    tuple val(meta), path("${meta.pair_id}.strelka2.unfiltered.snvs.vcf.gz"),   path("${meta.pair_id}.strelka2.unfiltered.snvs.vcf.gz.tbi"),   emit: vcf
    tuple val(meta), path("${meta.pair_id}.strelka2.unfiltered.indels.vcf.gz"), path("${meta.pair_id}.strelka2.unfiltered.indels.vcf.gz.tbi"), emit: indels

    script:
    def exome_flag        = params.genome ? '' : '--exome'
    def call_regions_flag = intervals_gz ? "--callRegions \"${intervals_gz}\"" : ''
    """
    configureStrelkaSomaticWorkflow.py \\
        --normalBam "${normal_bam}" \\
        --tumorBam  "${tumor_bam}" \\
        --referenceFasta "${fasta}" \\
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
