#!/bin/bash

set -euo pipefail

DOWNLOAD_DIR="/sc-projects/sc-proj-cc12-ag-romagnani/ref/imprint/gnomAD_download"
OUTDIR="/sc-projects/sc-proj-cc12-ag-romagnani/ref/imprint/gnomAD"
REF_FASTA="/sc-projects/sc-proj-cc12-ag-romagnani/ref/imprint/genome/GRCh38_no_alt_analysis_set.fasta"
CHRS=(chr{1..22} chrX chrY)

EXOMES_OUT="${OUTDIR}/gnomad.exomes.v4.1.1.af_only.vcf.gz"
GENOMES_OUT="${OUTDIR}/gnomad.genomes.v4.1.1.af_only.vcf.gz"
GERMLINE_OUT="${OUTDIR}/gnomad.germline_resource.v4.1.1.hg38.vcf.gz"
PILEUP_OUT="${OUTDIR}/gnomad.pileup_summaries.v4.1.1.hg38.vcf.gz"

mkdir -p "${DOWNLOAD_DIR}" "${OUTDIR}"
cd "${OUTDIR}"

source /opt/miniforge/etc/profile.d/conda.sh
conda activate gatk

gsutil -m cp -n gs://gcp-public-data--gnomad/release/4.1.1/vcf/* "${DOWNLOAD_DIR}"

# PASS-only, left-aligned, AC+AF+AN-stripped per-chr VCF
vcf_complete () {
    local vcf=$1
    [[ -s "${vcf}" && -s "${vcf}.tbi" ]] &&
        bcftools index -n "${vcf}" >/dev/null 2>&1
}

process_chr () {
    local dataset=$1 chr=$2
    local part="${OUTDIR}/gnomad.${dataset}.v4.1.1.af_only.${chr}.vcf.bgz"
    if vcf_complete "${part}"; then
        return
    fi
    rm -f -- "${part}" "${part}.tbi"

    local src="${DOWNLOAD_DIR}/${dataset}/gnomad.${dataset}.v4.1.1.sites.${chr}.vcf.bgz"
    local tmp="${part}.tmp.$$"
    [[ -s "${src}" ]] || { echo "Missing source VCF: ${src}" >&2; return 1; }
    bcftools view -f PASS "${src}" \
        | bcftools norm -f "${REF_FASTA}" \
        | bcftools annotate --remove ^INFO/AC,INFO/AF,INFO/AN -O z -o "${tmp}"
    bcftools index -t "${tmp}"
    mv "${tmp}" "${part}"
    mv "${tmp}.tbi" "${part}.tbi"
}
export -f process_chr vcf_complete
export DOWNLOAD_DIR OUTDIR REF_FASTA

build_track () {
    local dataset=$1 out=$2

    if vcf_complete "${out}"; then
        echo "[$(date)] ${out} already exists, skipping"
        return
    fi
    rm -f -- "${out}" "${out}.tbi"

    echo "[$(date)] Building ${out} from ${dataset}..."
    printf '%s\n' "${CHRS[@]}" | xargs -P "$(nproc)" -I{} bash -c 'process_chr "$1" "$2"' _ "${dataset}" {}

    local parts=()
    for chr in "${CHRS[@]}"; do
        parts+=("${OUTDIR}/gnomad.${dataset}.v4.1.1.af_only.${chr}.vcf.bgz")
    done

    local tmp="${out}.tmp.$$"
    bcftools concat -O z -o "${tmp}" "${parts[@]}"
    bcftools index -t "${tmp}"
    mv "${tmp}" "${out}"
    mv "${tmp}.tbi" "${out}.tbi"
}

build_track exomes  "${EXOMES_OUT}"
build_track genomes "${GENOMES_OUT}"

# Mutect2 uses the exome AF track as its germline population resource.
if ! vcf_complete "${GERMLINE_OUT}"; then
    rm -f -- "${GERMLINE_OUT}" "${GERMLINE_OUT}.tbi"
    echo "[$(date)] Linking ${GERMLINE_OUT} -> ${EXOMES_OUT}"
    ln -f "${EXOMES_OUT}" "${GERMLINE_OUT}"
    ln -f "${EXOMES_OUT}.tbi" "${GERMLINE_OUT}.tbi"
else
    echo "[$(date)] ${GERMLINE_OUT} already exists, skipping"
fi

# Pileup summaries use biallelic common SNPs (AF 0.01-0.2).
if ! vcf_complete "${PILEUP_OUT}"; then
    rm -f -- "${PILEUP_OUT}" "${PILEUP_OUT}.tbi"
    echo "[$(date)] Building ${PILEUP_OUT}..."
    pileup_tmp="${PILEUP_OUT}.tmp.$$"
    bcftools view -v snps -m2 -M2 -i 'INFO/AF>=0.01 && INFO/AF<=0.2' \
        -O z -o "${pileup_tmp}" "${EXOMES_OUT}"
    bcftools index -t "${pileup_tmp}"
    mv "${pileup_tmp}" "${PILEUP_OUT}"
    mv "${pileup_tmp}.tbi" "${PILEUP_OUT}.tbi"
else
    echo "[$(date)] ${PILEUP_OUT} already exists, skipping"
fi

rm -f "${OUTDIR}"/gnomad.exomes.v4.1.1.af_only.chr*.vcf.bgz*
rm -f "${OUTDIR}"/gnomad.genomes.v4.1.1.af_only.chr*.vcf.bgz*
