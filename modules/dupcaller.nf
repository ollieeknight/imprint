process DUPCALLER_TRIM_LANE {
    label 'process_low'
    label 'process_long'
    tag "${meta.id}"
    container "${params.container_dupcaller}"
    publishDir { "${params.outdir}/${meta.donor}/samples/${meta.cell_type}/qc/dupcaller" }, mode: 'copy', pattern: '*_barcode_metrics.json'

    input:
    tuple val(meta), path(reads1), path(reads2)
    path(allowlist)
    path(trim_script)

    output:
    tuple val(meta), path("${meta.id}_${meta.run_id}_R1_dupcaller.fastq.gz"), path("${meta.id}_${meta.run_id}_R2_dupcaller.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}_${meta.run_id}_barcode_metrics.json"), emit: metrics

    script:
    """
    python "${trim_script}" \
        --read1 "${reads1}" \
        --read2 "${reads2}" \
        --output1 "${meta.id}_${meta.run_id}_R1_dupcaller.fastq.gz" \
        --output2 "${meta.id}_${meta.run_id}_R2_dupcaller.fastq.gz" \
        --allowlist "${allowlist}" \
        --metrics "${meta.id}_${meta.run_id}_barcode_metrics.json"
    """
}

process DUPCALLER_ALIGN_LANE {
    label 'process_dynamic'
    label 'process_very_long'
    tag "${meta.id}"
    container "${params.container_align}"
    memory { 48.GB }
    cpus   { (reads1.size() + reads2.size()) < 5.GB ? 8 : (reads1.size() + reads2.size()) < 15.GB ? 12 : (reads1.size() + reads2.size()) < 60.GB ? 16 : 24 }

    input:
    tuple val(meta), path(reads1), path(reads2)

    output:
    tuple val(meta), path("${meta.id}_${meta.run_id}.dupcaller.unsorted.bam"), emit: bam

    script:
    """
    bwa-mem3 mem -C -K 100000000 -t ${task.cpus} \
        -R "@RG\\tID:${meta.id}.${meta.run_id}\\tSM:${meta.id}\\tPL:ILLUMINA\\tLB:${meta.id}" \
        "${params.bwa_mem3_index}" "${reads1}" "${reads2}" \
      | samtools view -1 -@ ${task.cpus} -o "${meta.id}_${meta.run_id}.dupcaller.unsorted.bam"
    """
}

process DUPCALLER_SORT_LANE {
    label 'process_dynamic'
    label 'process_very_long'
    tag "${meta.id}"
    container "${params.container_samtools}"
    memory { bam.size() < 20.GB ? 24.GB : 48.GB }
    cpus   { 8 }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${bam.name.replace('.unsorted.bam', '.bam')}"), emit: bam

    script:
    def lane_bam    = bam.name.replace('.unsorted.bam', '.bam')
    def total_gb    = task.memory ? task.memory.toGiga() : 24
    def sort_mem_gb = Math.max(1, ((total_gb - 4) / task.cpus).intValue())
    """
    samtools sort -@ ${task.cpus} -m ${sort_mem_gb}G -O bam -o "${lane_bam}" ${bam}
    """
}

process DUPCALLER_MARK_DUPLICATES {
    label 'process_high'
    label 'process_low_memory'
    label 'process_long'
    tag "${meta.id}"
    container "${params.container_gatk}"
    publishDir { "${params.outdir}/${meta.donor}/samples/${meta.cell_type}/qc/dupcaller" }, mode: 'copy', pattern: '*_dupcaller_markdup_metrics.txt'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.dupcaller.bam"), path("${meta.id}.dupcaller.bam.bai"), emit: bam
    tuple val(meta), path("${meta.id}_dupcaller_markdup_metrics.txt"), emit: metrics

    script:
    def mem_mb = task.memory ? (task.memory.mega * 0.8).intValue() : 12288
    """
    gatk --java-options "-Xmx${mem_mb}M" MarkDuplicates \
        --INPUT "${bam}" \
        --OUTPUT "${meta.id}.dupcaller.bam" \
        --METRICS_FILE "${meta.id}_dupcaller_markdup_metrics.txt" \
        --READ_NAME_REGEX '(?:.*:)?([0-9]+)[^:]*:([0-9]+)[^:]*:([0-9]+)[^:]*\$' \
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE ${params.optical_dup_dist} \
        --DUPLEX_UMI \
        --TAGGING_POLICY OpticalOnly \
        --BARCODE_TAG DB \
        --CREATE_INDEX true \
        --VALIDATION_STRINGENCY SILENT
    mv "${meta.id}.dupcaller.bai" "${meta.id}.dupcaller.bam.bai"
    """
}

process DUPCALLER_VALIDATE_BAM {
    label 'process_medium'
    label 'process_long'
    tag "${meta.id}"
    container "${params.container_align}"
    publishDir { "${params.outdir}/${meta.donor}/samples/${meta.cell_type}/qc/dupcaller" }, mode: 'copy', pattern: '*_dupcaller_tag_validation.tsv'

    input:
    tuple val(meta), path(bam), path(bai)
    path(allowlist)
    path(validator)

    output:
    tuple val(meta), path(bam), path(bai), emit: bam
    tuple val(meta), path("${meta.id}_dupcaller_tag_validation.tsv"), emit: metrics

    script:
    """
    samtools quickcheck -v "${bam}"
    samtools view -@ ${task.cpus} "${bam}" \
      | awk -v allowlist="${allowlist}" -f "${validator}" \
      > "${meta.id}_dupcaller_tag_validation.tsv"
    """
}

process DUPCALLER_VALIDATE_CRAM {
    label 'process_medium'
    label 'process_long'
    tag "${meta.id}"
    container "${params.container_align}"
    publishDir { "${params.outdir}/${meta.donor}/samples/${meta.cell_type}/qc/dupcaller" }, mode: 'copy'

    input:
    tuple val(meta), path(cram), path(crai)
    path(allowlist)
    path(validator)

    output:
    tuple val(meta), path("${meta.id}_dupcaller_cram_tag_validation.tsv"), emit: metrics

    script:
    """
    samtools quickcheck -v "${cram}"
    samtools view -@ ${task.cpus} -T "${params.ref_fasta}" "${cram}" \
      | awk -v allowlist="${allowlist}" -f "${validator}" \
      > "${meta.id}_dupcaller_cram_tag_validation.tsv"
    """
}

process DUPCALLER_CALL {
    label 'process_dynamic'
    label 'process_very_long'
    tag "${meta.pair_id}"
    container "${params.container_dupcaller}"
    publishDir { "${params.outdir}/${meta.donor}/pairs/${meta.pair_dir}/dupcaller" }, mode: 'copy'

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path(reference)
    path(reference_fai)
    path(reference_h5)
    path(trinuc_h5)
    path(homopolymer_h5)
    path(str_h5)
    path(dbs_h5)
    path(germline_vcf)
    path(germline_tbi)
    path(noise_masks)
    path(noise_mask_indexes)
    path(indel_epon)
    path(indel_epon_tbi)
    path(region_file)
    val(use_noise_masks)
    val(max_zero_qual_fraction)

    output:
    tuple val(meta), path('calls'), emit: calls

    script:
    def noiseArg  = use_noise_masks ? "-m ${noise_masks.join(' ')}" : ''
    def indelArg  = params.dupcaller_indel_epon ? "-id ${indel_epon}" : ''
    def rescueArg = params.dupcaller_rescue ? '--rescue' : ''
    def seedArg   = params.dupcaller_seed != null ? "--seed ${params.dupcaller_seed}" : ''
    """
    bgzip -c "${region_file}" > targets.bed.gz
    tabix -p bed targets.bed.gz

    DupCaller.py call \
        -b "${tumor_bam}" \
        -n "${normal_bam}" \
        -f "${reference}" \
        -g "${germline_vcf}" \
        -R targets.bed.gz \
        -o "${meta.pair_id}" \
        -p ${task.cpus} \
        -r ${params.dupcaller_regions} \
        ${noiseArg} ${indelArg} ${rescueArg} ${seedArg} \
        -gaf ${params.dupcaller_germline_af_cutoff} \
        -maf ${params.dupcaller_max_af} \
        -d ${params.dupcaller_min_normal_depth} \
        -z ${max_zero_qual_fraction} \
        -tt ${params.dupcaller_trim_template} \
        -tr ${params.dupcaller_trim_read} \
        -mq ${params.dupcaller_mapq} \
        -w ${params.dupcaller_window_size}

    rm -rf "${meta.pair_id}/tmp"

    for required in \
        "${meta.pair_id}_stats.txt" \
        "${meta.pair_id}_call_params.log" \
        "${meta.pair_id}_coverage.bed.gz" \
        "${meta.pair_id}_coverage.bed.gz.tbi" \
        "${meta.pair_id}_duplex_family_strand_composition.txt" \
        "SBS/${meta.pair_id}_sbs.vcf" \
        "SBS/${meta.pair_id}_trinuc_by_duplex_group.txt" \
        "INDEL/${meta.pair_id}_indel.vcf" \
        "INDEL/${meta.pair_id}_indel_by_duplex_group.txt" \
        "DBS/${meta.pair_id}_dbs.vcf" \
        "DBS/${meta.pair_id}_dbs_by_duplex_group.txt" \
        "ERROR/${meta.pair_id}.amp.tn.txt" \
        "ERROR/${meta.pair_id}.amp.hp.txt" \
        "ERROR/${meta.pair_id}.amp.str.txt" \
        "ERROR/${meta.pair_id}.dmg.tn.txt" \
        "ERROR/${meta.pair_id}.dmg.hp.txt" \
        "ERROR/${meta.pair_id}.dmg.str.txt"; do
        test -s "${meta.pair_id}/\$required"
    done

    for vcf in \
        "${meta.pair_id}/SBS/${meta.pair_id}_sbs.vcf" \
        "${meta.pair_id}/INDEL/${meta.pair_id}_indel.vcf" \
        "${meta.pair_id}/DBS/${meta.pair_id}_dbs.vcf"; do
        awk 'BEGIN { i = 0 }
             /^##contig=<ID=/ { match(\$0, /ID=[^,>]+/); ord[substr(\$0, RSTART + 3, RLENGTH - 3)] = i++ }
             /^#/ { print "-1\\t" \$0; next }
             { print ord[\$1] "\\t" \$0 }' "\$vcf" \
          | sort -s -k1,1n -k3,3n \
          | cut -f2- \
          | bgzip -c > "\$vcf.gz"
        tabix -p vcf "\$vcf.gz"
    done

    mv "${meta.pair_id}" calls
    """
}

process DUPCALLER_ESTIMATE {
    label 'process_medium'
    label 'process_long'
    tag "${meta.pair_id}"
    container "${params.container_dupcaller}"
    publishDir { "${params.outdir}/${meta.donor}/pairs/${meta.pair_dir}/dupcaller" }, mode: 'copy',
        saveAs: { fn -> fn == 'burden' ? 'burden' : null }

    input:
    tuple val(meta), path(calls)
    path(reference)
    path(reference_h5)
    path(trinuc_h5)
    path(homopolymer_h5)
    path(str_h5)
    path(dbs_h5)

    output:
    tuple val(meta), path('burden'),            emit: burden
    tuple val(meta), path("${meta.pair_id}"),   emit: sample_dir

    script:
    """
    # sigProfilerPlotting caches its matplotlib figure templates inside its own
    # package directory, which is read-only in the image. Its own environment
    # variable redirects that cache; a task-local path also keeps concurrent
    # pairs from writing the same pickle.
    export SIGPROFILERPLOTTING_VOLUME="\$PWD/.sigprofilerplotting"

    cp -rL "${calls}" "${meta.pair_id}"
    find "${meta.pair_id}" -type f | sort > before.txt

    DupCaller.py estimate \
        -i "${meta.pair_id}" \
        -f "${reference}" \
        -r ${params.dupcaller_regions}

    find "${meta.pair_id}" -type f | sort > after.txt
    comm -13 before.txt after.txt > new.txt
    test -s new.txt

    mkdir -p burden
    while IFS= read -r produced; do
        relative="\${produced#${meta.pair_id}/}"
        mkdir -p "burden/\$(dirname "\$relative")"
        cp "\$produced" "burden/\$relative"
    done < new.txt

    for required in \
        "SBS/${meta.pair_id}_sbs_burden.txt" \
        "SBS/${meta.pair_id}_sbs_96_corrected.txt" \
        "INDEL/${meta.pair_id}_indel_burden.txt" \
        "DBS/${meta.pair_id}_dbs_burden.txt" \
        "${meta.pair_id}_duplex_allele_counts.txt"; do
        test -s "burden/\$required"
    done
    """
}

process DUPCALLER_SUMMARIZE {
    label 'process_low'
    tag 'cohort'
    container "${params.container_dupcaller}"
    publishDir "${params.outdir}/cohort/dupcaller", mode: 'copy'

    input:
    path(sample_dirs)

    output:
    path('dupcaller_cohort_summary.txt'), emit: summary
    path('dupcaller_cohort_summary_SBS96_*.txt'), emit: sbs96

    script:
    def directories = (sample_dirs instanceof List ? sample_dirs : [sample_dirs])
        .collect { dir -> "\"${dir.name}\"" }.join(' ')
    """
    DupCaller.py summarize -i ${directories} -o dupcaller_cohort_summary.txt
    """
}
