process VARLOCIRAPTOR_ALIGNMENT_PROPERTIES {
    label 'process_low'
    label 'process_low_memory'
    label 'process_very_short'
    tag "${sample_id}"
    container "${params.container_varlociraptor}"

    input:
    tuple val(meta), val(role), val(sample_id), path(bam), path(bai)

    output:
    tuple val(meta), val(role), val(sample_id), path("${sample_id}.alignment_properties.json"), emit: json

    script:
    """
    varlociraptor estimate alignment-properties "${params.ref_fasta}" --bams "${bam}" > "${sample_id}.alignment_properties.json"
    """
}

process VARLOCIRAPTOR_PREPROCESS {
    label 'process_medium'
    label 'process_low_memory'
    tag "${sample_id}"
    container "${params.container_varlociraptor}"

    input:
    tuple val(meta), val(role), val(sample_id), path(bam), path(bai), path(candidates_vcf), path(candidates_tbi), path(alignment_properties_json)

    output:
    tuple val(meta), val(role), path("${sample_id}.obs.bcf"), emit: bcf

    script:
    """
    varlociraptor preprocess variants "${params.ref_fasta}" \\
        --bam "${bam}" \\
        --candidates "${candidates_vcf}" \\
        --output "${sample_id}.obs.bcf" \\
        --alignment-properties "${alignment_properties_json}" \\
        --pairhmm-mode exact \\
        --atomic-candidate-variants \\
        --max-depth 500
    """
}

process VARLOCIRAPTOR_CALLVARIANTS {
    label 'process_medium'
    label 'process_low_memory'
    tag "${meta.pair_id}"
    container "${params.container_varlociraptor}"

    input:
    tuple val(meta), path(tumour_obs_bcf), path(normal_obs_bcf), path(scenario_yaml)

    output:
    tuple val(meta), path("${meta.pair_id}.varlociraptor.bcf"), emit: bcf

    script:
    """
    varlociraptor call variants --output "${meta.pair_id}.varlociraptor.bcf" generic \\
        --scenario "${scenario_yaml}" \\
        --obs tumor="${tumour_obs_bcf}" normal="${normal_obs_bcf}"
    """
}

process VARLOCIRAPTOR_INDEX {
    label 'process_low'
    label 'process_low_memory'
    label 'process_very_short'
    tag "${meta.pair_id}"
    container "${params.container_bcftools}"
    publishDir { "${params.outdir}/${meta.donor}/pairs/${meta.pair_dir}/variant_calling/raw" }, mode: 'copy'

    input:
    tuple val(meta), path(varlociraptor_bcf)

    output:
    tuple val(meta), path("${meta.pair_id}.varlociraptor.bcf"), path("${meta.pair_id}.varlociraptor.bcf.csi"), emit: bcf

    script:
    """
    bcftools index "${varlociraptor_bcf}"
    """
}

process VARLOCIRAPTOR_MERGE {
    label 'process_low'
    label 'process_low_memory'
    label 'process_very_short'
    tag "${meta.pair_id}"
    container "${params.container_bcftools}"
    publishDir { "${params.outdir}/${meta.donor}/pairs/${meta.pair_dir}/variant_calling" }, mode: 'copy'

    input:
    tuple val(meta), path(vep_vcf), path(vep_tbi), path(varlociraptor_bcf), path(varlociraptor_csi)

    output:
    tuple val(meta), path("${meta.pair_id}.somatic.vcf.gz"), path("${meta.pair_id}.somatic.vcf.gz.tbi"), emit: vcf

    script:
    // Varlociraptor represents deletions >50 bp as <DEL>; recover their explicit
    // alleles from the unchanged candidate/VEP VCF using CHROM, POS and SVLEN.
    """
    fields='PROB_SOMATIC_TUMOR_LOW,PROB_SOMATIC_TUMOR_HIGH,PROB_SHARED_CLONAL,PROB_SOMATIC_NORMAL_ONLY,PROB_GERMLINE_HET,PROB_GERMLINE_HET_GAIN,PROB_GERMLINE_HOM,PROB_ABSENT,PROB_ARTIFACT'
    columns=\$(printf 'INFO/%s,' \${fields//,/ })
    columns=\${columns%,}

    bcftools annotate -a "${varlociraptor_bcf}" -c "INFO/HINTS,\$columns" \\
        "${vep_vcf}" -O z -o exact.vcf.gz
    tabix -p vcf exact.vcf.gz

    bcftools query -i 'ALT="<DEL>"' \\
        -f "%CHROM\\t%POS\\t%INFO/SVLEN\\t%INFO/HINTS\\t%INFO/\${fields//,/\\\\t%INFO/}\\n" \\
        "${varlociraptor_bcf}" > symbolic_deletions.tsv
    bcftools query -i 'strlen(REF)>strlen(ALT)' -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' \\
        "${vep_vcf}" > explicit_deletions.tsv

    awk 'BEGIN { FS=OFS="\\t" }
        NR==FNR { key=\$1 FS \$2 FS (-\$3); if (key in prob) exit 42; prob[key]=\$4; for (i=5; i<=NF; i++) prob[key]=prob[key] OFS \$i; next }
        { key=\$1 FS \$2 FS (length(\$3)-length(\$4)); if (key in prob) print \$0,prob[key] }' \\
        symbolic_deletions.tsv explicit_deletions.tsv > deletion_probabilities.tsv
    bgzip deletion_probabilities.tsv
    tabix -s 1 -b 2 -e 2 deletion_probabilities.tsv.gz

    bcftools annotate -a deletion_probabilities.tsv.gz \\
        -c "CHROM,POS,REF,ALT,INFO/HINTS,\$columns" exact.vcf.gz \\
        -O z -o "${meta.pair_id}.somatic.vcf.gz"
    tabix -p vcf "${meta.pair_id}.somatic.vcf.gz"

    missing=\$(bcftools query -f "%INFO/HINTS\\t%INFO/\${fields//,/\\\\t%INFO/}\\n" "${meta.pair_id}.somatic.vcf.gz" \\
        | awk '{ for (i=2; i<=NF; i++) if (\$i==".") { if (\$1 !~ /(^|,)missing-data(,|\$)/) n++; break } } END { print n+0 }')
    test "\$missing" -eq 0 || { echo "ERROR: \$missing records lack VarLociraptor posteriors" >&2; exit 1; }
    """
}
