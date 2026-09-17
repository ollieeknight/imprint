process MERGE_DONOR_BAMS {
    label 'process_medium'
    label 'process_low_memory'
    label 'process_long'
    tag "${meta.donor}"
    container "${params.container_samtools}"

    input:
    tuple val(meta), path(bams), path(bais)

    output:
    tuple val(meta), path("${meta.donor}.all.merged.bam"), path("${meta.donor}.all.merged.bam.bai"), emit: bam

    script:
    """
    samtools merge -f -@ ${task.cpus} ${meta.donor}.all.merged.bam ${bams}
    samtools index ${meta.donor}.all.merged.bam
    """
}

process OPTITYPE_EXTRACT {
    label 'process_high'
    label 'process_low_memory'
    label 'process_short'
    tag "${meta.donor}"
    container "${params.container_yara}"

    input:
    tuple val(meta), path(bam), path(bai)
    path hla_files

    output:
    tuple val(meta), path("r1.mapped.bam"), path("r1.mapped.bam.bai"),
                     path("r2.mapped.bam"), path("r2.mapped.bam.bai"), emit: bams

    script:
    """
    # The donor BAM is duplicate-marked, not duplicate-removed. Exclude flagged
    # PCR/optical duplicates so HLA support reflects independent templates.
    samtools view -@ ${task.cpus} -h -F 0xD00 -b "${bam}" | \\
        samtools collate -@ ${task.cpus} -O - . | \\
        samtools fastq -@ ${task.cpus} \\
            -1 r1.fq.gz -2 r2.fq.gz \\
            -0 /dev/null -s /dev/null -N

    yara_mapper -t ${task.cpus} -f bam hla_reference_dna.fasta r1.fq.gz r2.fq.gz > mapped.bam

    samtools view -@ ${task.cpus} -hF 4 -f 0x40 -b mapped.bam | samtools sort -@ ${task.cpus} > r1.mapped.bam
    samtools view -@ ${task.cpus} -hF 4 -f 0x80 -b mapped.bam | samtools sort -@ ${task.cpus} > r2.mapped.bam
    samtools index r1.mapped.bam
    samtools index r2.mapped.bam
    """
}

process OPTITYPE_GENOTYPE {
    label 'process_medium'
    tag "${meta.donor}"
    container "${params.container_optitype}"
    publishDir { "${params.outdir}/${meta.donor}/hla" }, mode: 'copy'

    input:
    tuple val(meta), path(r1_bam), path(r1_bai), path(r2_bam), path(r2_bai)

    output:
    tuple val(meta), path("${meta.donor}_result.tsv"),        emit: result
    tuple val(meta), path("${meta.donor}_coverage_plot.pdf"), emit: plot, optional: true

    script:
    """
    printf '[mapping]\\nrazers3=razers3\\nthreads=${task.cpus}\\n[ilp]\\nsolver=glpk\\nthreads=1\\n[behavior]\\ndeletebam=true\\nunpaired_weight=0\\nuse_discordant=false\\n' > config.ini

    # Unset proxy variables to prevent Python's configparser DuplicateOptionError
    unset HTTP_PROXY HTTPS_PROXY http_proxy https_proxy

    OptiTypePipeline.py \\
        -i ${r1_bam} ${r2_bam} \\
        -c config.ini \\
        --dna \\
        --prefix "${meta.donor}" \\
        --outdir .
    """
}

process KIR_MAPPER {
    label 'process_high'
    label 'process_high_memory'
    label 'process_very_long'
    tag "${meta.donor}"
    container "${params.container_kir_mapper}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("ncopy"), path("genotype"), emit: raw

    script:
    def exome_flag = params.genome ? '' : '--exome'
    """
    # Write kir-mapper config with resolved db path; set HOME so kir-mapper finds it.
    cp /opt/.kir-mapper .kir-mapper
    sed -i "s|__KIR_MAPPER_DB__|${params.kirmapper_db}|g" .kir-mapper
    export HOME=\$PWD

    # Subset to KIR/LILRB locus + HLA-E/G locus + all alt/unplaced contigs + unmapped.
    # Regions:
    #   chr19:54000000-56000000  KIR + LILRB genes (LRC)
    #   chr6:29600000-30000000   HLA-E and HLA-G
    #   ALL_ALTS                 every non-primary contig (alt, unplaced, decoy); B-haplotype
    #                            KIR reads absent from GRCh38 primary may land anywhere here.
    # Two-pass for mate rescue: reads whose mate maps to these regions but the read itself
    # landed elsewhere are captured by name-based second pass.
    # -f 12 = both reads unmapped; avoids double-counting mates already rescued by name.
    ALL_ALTS=\$(samtools view -H "${bam}" | \
      awk '/^@SQ/ {sub(/.*SN:/, "", \$2); print \$2}' | \
      grep -vE '^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)\$' | \
      tr '\\n' ' ')

    samtools view -F 0xD00 "${bam}" \
      chr19:54000000-56000000 chr6:29600000-30000000 \
      \${ALL_ALTS} | \
      awk '{print \$1}' | sort -u > kir_qnames.txt

    { samtools view -H "${bam}" ; \
      samtools view -F 0xD00 -N kir_qnames.txt "${bam}" ; \
      samtools view -f 12 -F 0xD00 "${bam}" ; } \
      | samtools sort -@ ${task.cpus} -o kir_sorted.bam
    samtools index kir_sorted.bam

    kir-mapper map \\
        -bam kir_sorted.bam \\
        -db ${params.kirmapper_db} \\
        -threads ${task.cpus} \\
        ${exome_flag} \\
        -output results/

    kir-mapper ncopy \\
        -db ${params.kirmapper_db} \\
        -output results/ \\
        -threads ${task.cpus} \\
        ${exome_flag}

    kir-mapper genotype \\
        -db ${params.kirmapper_db} \\
        -output results/ \\
        -threads ${task.cpus} \\
        ${exome_flag}

    # Clean up intermediate SAM files from kir-mapper processing
    find results/ -name "*.sam*" -type f -delete
    mv results/ncopy .
    if [ -d results/genotype ]; then
        mv results/genotype .
    else
        mkdir -p genotype
    fi
    """
}

process KIR_COLLATE {
    label 'process_single'
    tag "${meta.donor}"
    container "${params.container_kir_mapper}"
    publishDir { "${params.outdir}/${meta.donor}/kir" }, mode: 'copy'

    input:
    tuple val(meta), path(ncopy), path(genotype)
    path collator

    output:
    tuple val(meta), path("${meta.donor}.kir.copy_number.tsv"),               emit: copy_number
    tuple val(meta), path("${meta.donor}.kir.calls.tsv"),                     emit: calls
    tuple val(meta), path("${meta.donor}.kir.genotype_candidates.tsv.gz"),    emit: candidates
    tuple val(meta), path("${meta.donor}.kir-mapper.raw.tar.gz"),             emit: raw_archive

    script:
    """
    python3 "${collator}" \\
        --donor "${meta.donor}" \\
        --ncopy-dir "${ncopy}" \\
        --genotype-dir "${genotype}" \\
        --output-dir .
    """
}

process PATHSEQ {
    label 'process_high'
    label 'process_high_cpu'
    label 'process_long'
    tag "${meta.id}"
    container "${params.container_gatk}"
    publishDir {
        "${params.outdir}/${meta.donor}/samples/${meta.cell_type}/pathseq"
    }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    tuple path(host_img), path(host_hss), path(microbe_img), path(microbe_dict), path(taxonomy_db)

    output:
    tuple val(meta), path("${meta.id}.pathseq.bam"),                  emit: bam
    tuple val(meta), path("${meta.id}.pathseq.scores.txt"),           emit: scores
    tuple val(meta), path("${meta.id}.pathseq.filter.metrics.txt"),   emit: filter_metrics
    tuple val(meta), path("${meta.id}.pathseq.score.metrics.txt"),    emit: score_metrics
    tuple val(meta), path("${meta.id}.pathseq.score.warnings.txt"),   emit: score_warnings

    script:
    def mem_gb       = task.memory ? task.memory.toGiga() - 8 : 8
    """
    # PathSeqFilterSpark clears the sequence dictionary before validating
    # --ignore-alignment-contigs with --is-host-aligned. Unmap EBV reads here
    # to avoid "Ignored sequence ... not found in input header".
    ebv_contig="\$(samtools view -H "${bam}" 2>/dev/null | awk -F '\\t' '
        \$1 == "@SQ" {
            for (i = 2; i <= NF; i++) {
                if (\$i == "SN:chrEBV" || \$i == "SN:EBV") {
                    sub(/^SN:/, "", \$i)
                    print \$i
                    exit 0
                }
            }
        }
    ')"

    if [ -n "\${ebv_contig}" ]; then
        echo "INFO: Unmapping reads aligned to \${ebv_contig} to bypass GATK bug..." >&2
        samtools view -h "${bam}" | awk -v contig="\${ebv_contig}" '
            BEGIN { OFS="\\t" }
            function testbit(f,b){return int(f/b)%2}
            function setbit(f,b){return testbit(f,b)?f:f+b}
            function clrbit(f,b){return testbit(f,b)?f-b:f}
            /^@/ { print; next }
            {
                if (\$3 == contig || \$7 == contig || (\$7 == "=" && \$3 == contig)) {
                    flag=\$2
                    flag=setbit(flag,4)                 # unmapped
                    if (testbit(flag,1)) flag=setbit(flag,8)  # mate unmapped (if paired)
                    flag=clrbit(flag,2); flag=clrbit(flag,16)
                    flag=clrbit(flag,32); flag=clrbit(flag,256); flag=clrbit(flag,2048)
                    \$2=flag; \$3="*"; \$4=0; \$5=0; \$6="*"; \$7="*"; \$8=0; \$9=0
                }
                print
            }
        ' | samtools view -b -o pathseq_input.bam -
    else
        echo "INFO: Neither chrEBV nor EBV is present in ${bam} header." >&2
        ln -s "${bam}" pathseq_input.bam
    fi

    gatk --java-options "-Xmx${mem_gb}g --add-opens java.base/javax.security.auth=ALL-UNNAMED" PathSeqPipelineSpark \
        --input pathseq_input.bam \
        --output "${meta.id}.pathseq.bam" \
        --filter-bwa-image "${host_img}" \
        --kmer-file "${host_hss}" \
        --microbe-bwa-image "${microbe_img}" \
        --microbe-dict "${microbe_dict}" \
        --taxonomy-file "${taxonomy_db}" \
        --scores-output "${meta.id}.pathseq.scores.txt" \
        --is-host-aligned true \
        --host-kmer-thresh 2 \
        --min-clipped-read-length 60 \
        --min-score-identity 0.90 \
        --identity-margin 0.02 \
        --filter-metrics "${meta.id}.pathseq.filter.metrics.txt" \
        --score-metrics "${meta.id}.pathseq.score.metrics.txt" \
        --score-warnings "${meta.id}.pathseq.score.warnings.txt" \
        --filter-duplicates true \
        --spark-master local[${task.cpus}]
    """
}

process TELSEQ {
    label 'process_low'
    label 'process_short'
    tag "${meta.id}"
    container "${params.container_telseq}"
    publishDir { "${params.outdir}/${meta.donor}/samples/${meta.cell_type}/telseq" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.telseq.txt"), optional: true, emit: telseq

    script:
    def exome_arg = (!params.genome && params.intervals_bed) ? "-e ${params.intervals_bed}" : ''
    def rlen      = meta.read_length ?: 100
    """
    telseq -m -r ${rlen} ${exome_arg} ${bam} > ${meta.id}.telseq.txt
    """
}
