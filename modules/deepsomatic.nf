process CALL {
    label 'process_high'
    label 'process_high_cpu'
    tag "${meta.pair_id}"
    container "${params.container_deepsomatic}"
    publishDir { "${params.outdir}/${meta.donor}/pairs/${meta.pair_dir}/variant_calling/raw" }, mode: 'copy'

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)

    output:
    tuple val(meta), path("${meta.pair_id}.deepsomatic.unfiltered.vcf.gz"), path("${meta.pair_id}.deepsomatic.unfiltered.vcf.gz.tbi"), emit: vcf

    script:
    def model_type    = params.genome ? 'WGS' : 'WES'
    def effective_bed = params.off_target ? params.padded_intervals_bed : params.intervals_bed
    def regions_flag  = params.genome ? '' : "--regions=\"${effective_bed}\""
    """
    mkdir -p intermediate logs

    run_deepsomatic \\
        --model_type=${model_type} \\
        --ref="${params.ref_fasta}" \\
        --reads_tumor="${tumor_bam}" \\
        --reads_normal="${normal_bam}" \\
        --output_vcf="${meta.pair_id}.deepsomatic.unfiltered.vcf.gz" \\
        --sample_name_tumor="${meta.tumor_id}" \\
        --sample_name_normal="${meta.normal_id}" \\
        ${regions_flag} \\
        --num_shards=${task.cpus} \\
        --logging_dir=logs \\
        --intermediate_results_dir=intermediate

    rm -rf intermediate logs
    """
}
