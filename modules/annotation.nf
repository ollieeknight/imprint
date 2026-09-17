process STRELKA2_MERGE {
    label 'process_low'
    label 'process_low_memory'
    label 'process_very_short'
    tag "${meta.pair_id}"
    container "${params.container_bcftools}"

    input:
    tuple val(meta), path(snv_vcf), path(snv_tbi), path(indel_vcf), path(indel_tbi)

    output:
    tuple val(meta), path("${meta.pair_id}.strelka2.all.vcf.gz"), path("${meta.pair_id}.strelka2.all.vcf.gz.tbi"), emit: all_vcf

    script:
    """
    bcftools concat -a "${snv_vcf}" "${indel_vcf}" \\
        | bcftools view -f 'PASS' \\
        | bcftools sort -O z -o "${meta.pair_id}.strelka2.all.vcf.gz"
    tabix -p vcf "${meta.pair_id}.strelka2.all.vcf.gz"
    """
}

process ENSEMBLE_CONSENSUS {
    label 'process_low'
    label 'process_low_memory'
    label 'process_very_short'
    tag "${meta.pair_id}"
    container "${params.container_bcftools}"

    input:
    tuple val(meta), path(mutect2_vcf), path(mutect2_tbi), path(strelka_vcf), path(strelka_tbi), path(deepsomatic_vcf), path(deepsomatic_tbi)

    output:
    tuple val(meta), path("${meta.pair_id}.consensus.vcf.gz"), path("${meta.pair_id}.consensus.vcf.gz.tbi"), emit: vcf

    script:
    """
    # Split multiallelics and atomise MNPs before comparing caller support.
    # Mutect2 emits phased MNPs; Strelka2 and DeepSomatic only emit single
    # substitutions, so an intact MNP would match neither. --old-rec-tag MNP
    # records the source record so the codon can be regrouped downstream.
    for spec in "m2:${mutect2_vcf}" "st:${strelka_vcf}" "ds:${deepsomatic_vcf}"; do
        prefix="\${spec%%:*}"
        source_vcf="\${spec#*:}"
        bcftools view -f 'PASS' --drop-genotypes "\$source_vcf" \\
            | bcftools norm -m-any -f "${params.ref_fasta}" \\
            | bcftools norm --atomize --old-rec-tag MNP -f "${params.ref_fasta}" \\
            | bcftools sort -O z -o "\${prefix}.sites.vcf.gz"
    done

    tabix -f -p vcf m2.sites.vcf.gz
    tabix -f -p vcf st.sites.vcf.gz
    tabix -f -p vcf ds.sites.vcf.gz

    # Deduplicate exact alleles while retaining distinct alleles at one position.
    bcftools concat -a m2.sites.vcf.gz st.sites.vcf.gz ds.sites.vcf.gz \
        | bcftools sort \
        | bcftools norm -d exact \
        | bcftools annotate -x INFO -O z -o union.sites.vcf.gz
    tabix -f -p vcf union.sites.vcf.gz

    cat <<EOF > caller_hdr.txt
##INFO=<ID=MUTECT2_CALL,Number=0,Type=Flag,Description="Site present in Mutect2 PASS output">
##INFO=<ID=STRELKA2_CALL,Number=0,Type=Flag,Description="Site present in Strelka2 PASS output">
##INFO=<ID=DEEPSOMATIC_CALL,Number=0,Type=Flag,Description="Site present in DeepSomatic PASS output">
EOF

    bcftools annotate -h caller_hdr.txt \
        -a m2.sites.vcf.gz -c CHROM,POS,REF,ALT,INFO/MNP --mark-sites +MUTECT2_CALL \
        union.sites.vcf.gz -O z -o union.m2.vcf.gz
    tabix -f -p vcf union.m2.vcf.gz

    bcftools annotate \\
        -a st.sites.vcf.gz -c CHROM,POS,REF,ALT --mark-sites +STRELKA2_CALL \\
        union.m2.vcf.gz -O z -o union.m2st.vcf.gz
    tabix -f -p vcf union.m2st.vcf.gz

    bcftools annotate \\
        -a ds.sites.vcf.gz -c CHROM,POS,REF,ALT --mark-sites +DEEPSOMATIC_CALL \\
        union.m2st.vcf.gz -O z -o union.m2stds.vcf.gz
    tabix -f -p vcf union.m2stds.vcf.gz

    bcftools view union.m2stds.vcf.gz | awk '
    BEGIN { OFS="\\t" }
    /^##/ { print; next }
    /^#CHROM/ {
        print "##INFO=<ID=CALLER_SUPPORT,Number=1,Type=String,Description=\\"Comma-separated list of variant callers supporting this site (alphabetical: deepsomatic,mutect2,strelka2)\\">"
        print; next
    }
    {
        support = ""
        if (\$8 ~ /(^|;)DEEPSOMATIC_CALL(;|\$)/)  support = "deepsomatic"
        if (\$8 ~ /(^|;)MUTECT2_CALL(;|\$)/)      support = (length(support) > 0 ? support "," : "") "mutect2"
        if (\$8 ~ /(^|;)STRELKA2_CALL(;|\$)/)     support = (length(support) > 0 ? support "," : "") "strelka2"
        \$8 = \$8 ";CALLER_SUPPORT=" support
        print
    }
    ' | bcftools view -O z -o "${meta.pair_id}.consensus.vcf.gz"

    tabix -f -p vcf "${meta.pair_id}.consensus.vcf.gz"
    """
}

process VAFATOR {
    label 'process_medium'
    label 'process_low_memory'
    tag "${meta.pair_id}"
    container "${params.container_vafator}"

    input:
    tuple val(meta), path(consensus_vcf), path(consensus_tbi), path(tumour_bam), path(tumour_bai), path(normal_bam), path(normal_bai)

    output:
    tuple val(meta), path("${meta.pair_id}.vaf.vcf"), emit: vcf

    script:
    def mq = 20
    def bq = 20
    """
    vafator \\
        --input-vcf "${consensus_vcf}" \\
        --output-vcf "${meta.pair_id}.vaf.vcf" \\
        --bam "${meta.tumor_id}" "${tumour_bam}" \\
        --bam "${meta.normal_id}" "${normal_bam}" \\
        --mapping-quality ${mq} \\
        --base-call-quality ${bq} \\
        --exclude-ambiguous-bases \\
        --num-processes ${task.cpus}
    """
}

process TUMOR_NORMAL_FILTER {
    label 'process_low'
    label 'process_low_memory'
    label 'process_very_short'
    tag "${meta.pair_id}"
    container "${params.container_bcftools}"


    input:
    tuple val(meta), path(vaf_vcf)

    output:
    tuple val(meta), path("${meta.pair_id}.somatic_filtered.vcf.gz"), path("${meta.pair_id}.somatic_filtered.vcf.gz.tbi"), emit: vcf

    script:
    """
    bgzip -c "${vaf_vcf}" > input.vcf.gz
    tabix -p vcf input.vcf.gz

    # Mark ON_TARGET sites if off-target mode is active (captures only non-padded intervals)
    if [ "${params.off_target}" = "true" ] && [ -n "${params.intervals_bed}" ]; then
        printf '##INFO=<ID=ON_TARGET,Number=0,Type=Flag,Description="Variant overlaps capture target intervals (non-padded)">\\n' > on_target_hdr.txt
        awk 'BEGIN{OFS="\\t"} !/^#/{print \$1, \$2+1, \$3}' "${params.intervals_bed}" | bgzip -c > on_target_annot.bed.gz
        tabix -s1 -b2 -e3 -c '#' on_target_annot.bed.gz
        bcftools annotate --mark-sites "+ON_TARGET" -a on_target_annot.bed.gz -c CHROM,FROM,TO -h on_target_hdr.txt input.vcf.gz -O z -o input.marked.vcf.gz
        tabix -p vcf input.marked.vcf.gz
        mv input.marked.vcf.gz input.vcf.gz
        tabix -f -p vcf input.vcf.gz
    fi

    # Pass through all PASS variants (ensemble filter; no statistical gating)
    bcftools filter -i "FILTER='PASS'" input.vcf.gz -O z -o "${meta.pair_id}.somatic_filtered.vcf.gz"
    tabix -p vcf "${meta.pair_id}.somatic_filtered.vcf.gz"
    """
}

process VEP_ANNOTATE {
    label 'process_high'
    label 'process_high_memory'
    tag "${meta.pair_id}"
    container "${params.container_vep}"
    publishDir { "${params.outdir}/${meta.donor}/pairs/${meta.pair_dir}/dupcaller/annotated" },
        mode: 'copy', enabled: params.dupcaller

    input:
    tuple val(meta), val(mutation_type), path(vcf), path(tbi)
    // Distinct staging names are required in SNV-only mode because the same
    // physical VCF/index satisfy both mandatory SpliceAI plugin arguments.
    path spliceai_snv_vcf,   stageAs: 'spliceai.snv.vcf.gz'
    path spliceai_snv_tbi,   stageAs: 'spliceai.snv.vcf.gz.tbi'
    path spliceai_indel_vcf, stageAs: 'spliceai.indel.vcf.gz'
    path spliceai_indel_tbi, stageAs: 'spliceai.indel.vcf.gz.tbi'
    path dbnsfp_gz,           stageAs: 'dbNSFP.gz'
    path dbnsfp_tbi,          stageAs: 'dbNSFP.gz.tbi'

    output:
    tuple val(meta), val(mutation_type), path("${meta.pair_id}.${mutation_type}.vcf.gz"), path("${meta.pair_id}.${mutation_type}.vcf.gz.tbi"), emit: vcf

    script:
    def cosmic_flag = params.cosmic_vcf ? "--custom ${params.cosmic_vcf},COSMIC,vcf,exact,0,ID" : ''
    def cosmic_noncoding_flag = params.cosmic_noncoding_vcf ? "--custom ${params.cosmic_noncoding_vcf},COSMIC_NONCODING,vcf,exact,0,ID" : ''
    def custom_flags = [
        "--custom ${params.gnomad_exomes_vep_vcf},gnomADv4e,vcf,exact,0,AF,AN",
        "--custom ${params.gnomad_genomes_vep_vcf},gnomADv4g,vcf,exact,0,AF,AN",
        cosmic_flag,
        cosmic_noncoding_flag,
    ].findAll { v -> v }.join(" \\\n        ")
    """
    if [ \$(zcat "${vcf}" | grep -v "^#" | wc -l) -eq 0 ]; then
        cp "${vcf}" "${meta.pair_id}.${mutation_type}.vcf.gz"
        tabix -f -p vcf "${meta.pair_id}.${mutation_type}.vcf.gz"
        exit 0
    fi

    vep \\
        -i "${vcf}" \\
        -o "${meta.pair_id}.vep.vcf.gz" \\
        --vcf --compress_output bgzip --pick \\
        --offline --dir_cache "${params.vep_cache}" \\
        --fasta "${params.ref_fasta}" \\
        --species homo_sapiens \\
        --assembly GRCh38 \\
        --sift b \\
        --polyphen b \\
        --hgvs \\
        --symbol \\
        --domains \\
        --regulatory \\
        --canonical \\
        --protein \\
        --biotype \\
        --variant_class \\
        --mane \\
        --dir_plugins "${params.vep_plugins_dir}" \\
        --plugin AlphaMissense,file="${params.alphamissense_tsv}" \\
        --plugin SpliceAI,snv="${spliceai_snv_vcf}",indel="${spliceai_indel_vcf}" \\
        --plugin dbNSFP,${dbnsfp_gz},REVEL_score \\
        --plugin NMD \\
        --plugin SpliceRegion \\
        --plugin pLI \\
        ${custom_flags} \\
        --fork ${task.cpus}

    mv "${meta.pair_id}.vep.vcf.gz" "${meta.pair_id}.${mutation_type}.vcf.gz"
    tabix -f -p vcf "${meta.pair_id}.${mutation_type}.vcf.gz"
    """
}
