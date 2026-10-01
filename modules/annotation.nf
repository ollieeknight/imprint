process STRELKA2_MERGE {
    label 'process_low'
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
    tag "${meta.pair_id}"
    container "${params.container_bcftools}"

    input:
    tuple val(meta), path(mutect2_vcf), path(mutect2_tbi), path(strelka_vcf), path(strelka_tbi), path(deepsomatic_vcf), path(deepsomatic_tbi)
    tuple path(fasta), path(fai), path(dict)

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
            | bcftools norm -m-any -f "${fasta}" \\
            | bcftools norm --atomize --old-rec-tag MNP -f "${fasta}" \\
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
    tag "${meta.pair_id}"
    container "${params.container_vafator}"

    input:
    tuple val(meta), path(consensus_vcf), path(consensus_tbi), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)

    output:
    tuple val(meta), path("${meta.pair_id}.vaf.vcf"), emit: vcf

    script:
    """
    vafator \\
        --input-vcf "${consensus_vcf}" \\
        --output-vcf "${meta.pair_id}.vaf.vcf" \\
        --bam "${meta.tumor_id}" "${tumor_bam}" \\
        --bam "${meta.normal_id}" "${normal_bam}" \\
        --mapping-quality 20 \\
        --base-call-quality 20 \\
        --exclude-ambiguous-bases \\
        --num-processes ${task.cpus}
    """
}

process MARK_ON_TARGET {
    label 'process_low'
    tag "${meta.pair_id}"
    container "${params.container_bcftools}"

    input:
    tuple val(meta), path(vaf_vcf)
    path(target_bed)

    output:
    tuple val(meta), path("${meta.pair_id}.vaf.vcf.gz"), path("${meta.pair_id}.vaf.vcf.gz.tbi"), emit: vcf

    script:
    // Off-target runs call on padded intervals; flag sites inside the unpadded targets.
    // Otherwise this only compresses and indexes Vafator's output for VEP.
    if (target_bed) {
        """
        printf '##INFO=<ID=ON_TARGET,Number=0,Type=Flag,Description="Variant overlaps capture target intervals (non-padded)">\\n' > on_target_hdr.txt
        awk 'BEGIN{OFS="\\t"} !/^#/{print \$1, \$2+1, \$3}' "${target_bed}" | bgzip -c > on_target.bed.gz
        tabix -s1 -b2 -e3 -c '#' on_target.bed.gz
        bcftools annotate --mark-sites +ON_TARGET -a on_target.bed.gz -c CHROM,FROM,TO -h on_target_hdr.txt \\
            "${vaf_vcf}" -O z -o "${meta.pair_id}.vaf.vcf.gz"
        tabix -p vcf "${meta.pair_id}.vaf.vcf.gz"
        """
    } else {
        """
        bgzip -c "${vaf_vcf}" > "${meta.pair_id}.vaf.vcf.gz"
        tabix -p vcf "${meta.pair_id}.vaf.vcf.gz"
        """
    }
}

process VEP_ANNOTATE {
    label 'process_high'
    tag "${meta.pair_id}"
    container "${params.container_vep}"
    publishDir { "${params.outdir}/${meta.donor}/pairs/${meta.pair_dir}/dupcaller/annotated" },
        mode: 'copy', enabled: params.dupcaller

    input:
    tuple val(meta), val(mutation_type), path(vcf), path(tbi)
    tuple path(fasta), path(fai), path(dict)
    tuple path(vep_cache), path(vep_plugins)
    // Fixed names: in SNV-only mode the SNV VCF fills both SpliceAI arguments.
    tuple path(spliceai_snv, stageAs: 'spliceai.snv.vcf.gz'), path(spliceai_snv_tbi, stageAs: 'spliceai.snv.vcf.gz.tbi'),
          path(spliceai_indel, stageAs: 'spliceai.indel.vcf.gz'), path(spliceai_indel_tbi, stageAs: 'spliceai.indel.vcf.gz.tbi')
    tuple path(dbnsfp, stageAs: 'dbNSFP.gz'), path(dbnsfp_tbi, stageAs: 'dbNSFP.gz.tbi')
    tuple path(alphamissense), path(alphamissense_tbi)
    tuple path(gnomad_exomes), path(gnomad_exomes_tbi), path(gnomad_genomes), path(gnomad_genomes_tbi)
    tuple path(cosmic), path(cosmic_tbi), path(cosmic_noncoding), path(cosmic_noncoding_tbi)

    output:
    tuple val(meta), val(mutation_type), path("${meta.pair_id}.${mutation_type}.vcf.gz"), path("${meta.pair_id}.${mutation_type}.vcf.gz.tbi"), emit: vcf

    script:
    def custom_flags = [
        "--custom ${gnomad_exomes},gnomADv4e,vcf,exact,0,AF,AN",
        "--custom ${gnomad_genomes},gnomADv4g,vcf,exact,0,AF,AN",
        cosmic ? "--custom ${cosmic},COSMIC,vcf,exact,0,ID" : '',
        cosmic_noncoding ? "--custom ${cosmic_noncoding},COSMIC_NONCODING,vcf,exact,0,ID" : '',
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
        --offline --dir_cache "${vep_cache}" \\
        --fasta "${fasta}" \\
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
        --dir_plugins "${vep_plugins}" \\
        --plugin AlphaMissense,file="${alphamissense}" \\
        --plugin SpliceAI,snv="${spliceai_snv}",indel="${spliceai_indel}" \\
        --plugin dbNSFP,${dbnsfp},REVEL_score \\
        --plugin NMD \\
        --plugin SpliceRegion \\
        --plugin pLI \\
        ${custom_flags} \\
        --fork ${task.cpus}

    mv "${meta.pair_id}.vep.vcf.gz" "${meta.pair_id}.${mutation_type}.vcf.gz"
    tabix -f -p vcf "${meta.pair_id}.${mutation_type}.vcf.gz"
    """
}
