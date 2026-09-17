#!/bin/bash
#SBATCH --job-name=pathseq_ref
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=/sc-projects/sc-proj-cc12-ag-romagnani/ref/imprint/characterisation/pathseq/pathseq_reference_%j.log
#SBATCH --error=/sc-projects/sc-proj-cc12-ag-romagnani/ref/imprint/characterisation/pathseq/pathseq_reference_%j.log

set -euo pipefail

HOST_FASTA="${HOST_FASTA:-/sc-projects/sc-proj-cc12-ag-romagnani/ref/imprint/genome/GRCh38_no_alt_analysis_set.fasta}"
OUTDIR="${OUTDIR:-/sc-projects/sc-proj-cc12-ag-romagnani/ref/imprint/characterisation/pathseq}"
GATK_VERSION_REQUIRED="${GATK_VERSION_REQUIRED:-4.6.2.0}"

HOST_NO_EBV="${OUTDIR}/host/GRCh38_no_alt_analysis_set.no_chrEBV.fasta"
HOST_NO_EBV_IMG="${HOST_NO_EBV}.img"
HOST_NO_EBV_HSS="${OUTDIR}/host/GRCh38_no_alt_analysis_set.no_chrEBV.hss"
MICROBE_FASTA="${OUTDIR}/microbe/viral_refseq.fasta"
MICROBE_DICT="${OUTDIR}/microbe/viral_refseq.dict"
MICROBE_IMG="${MICROBE_FASTA}.img"
TAXDUMP="${OUTDIR}/taxonomy/taxdump.current.tar.gz"
REFSEQ_CATALOG="${OUTDIR}/taxonomy/refseq.current.catalog.gz"
TAXONOMY_DB="${OUTDIR}/taxonomy/viral_refseq.pathseq_taxonomy.db"

log() { echo "[$(date --iso-8601=seconds)] $*"; }
die() { log "ERROR: $*" >&2; exit 1; }

[[ -f "${HOST_FASTA}" ]] || die "Primary host FASTA not found: ${HOST_FASTA}"

source /opt/miniforge/etc/profile.d/conda.sh
conda activate gatk

for command_name in datasets unzip samtools gatk wget; do
    command -v "${command_name}" >/dev/null 2>&1 || die "Required command not found: ${command_name}"
done

gatk_version="$(gatk --version 2>&1)"
grep -Fq "${GATK_VERSION_REQUIRED}" <<< "${gatk_version}" || \
    die "GATK ${GATK_VERSION_REQUIRED} is required to match pipeline runtime; found: ${gatk_version}"

mkdir -p "${OUTDIR}"/{host,microbe,taxonomy}

for output in \
    "${HOST_NO_EBV}" "${HOST_NO_EBV}.fai" "${HOST_NO_EBV_IMG}" "${HOST_NO_EBV_HSS}" \
    "${MICROBE_FASTA}" "${MICROBE_FASTA}.fai" "${MICROBE_DICT}" "${MICROBE_IMG}" \
    "${TAXDUMP}" "${REFSEQ_CATALOG}" "${TAXONOMY_DB}"; do
    [[ ! -e "${output}" ]] || die "Output already exists: ${output}"
done

build_tmp="$(mktemp -d "${SLURM_TMPDIR:-/tmp}/pathseq-reference.XXXXXX")"
trap 'rm -rf -- "${build_tmp}"' EXIT

log "Building PathSeq host reference from GRCh38 minus EBV"
[[ -f "${HOST_FASTA}.fai" ]] || samtools faidx "${HOST_FASTA}"
ebv_contig="$(awk '$1 == "chrEBV" || $1 == "EBV" { print $1 }' "${HOST_FASTA}.fai")"
[[ -n "${ebv_contig}" ]] || die "Neither chrEBV nor EBV contig found in ${HOST_FASTA}.fai"
log "Detected EBV contig in host FASTA: ${ebv_contig}"

awk -v ebv="${ebv_contig}" '$1 != ebv { print $1 }' "${HOST_FASTA}.fai" > "${build_tmp}/host_contigs.txt"
samtools faidx "${HOST_FASTA}" -r "${build_tmp}/host_contigs.txt" > "${HOST_NO_EBV}"
samtools faidx "${HOST_NO_EBV}"

original_contigs="$(wc -l < "${HOST_FASTA}.fai")"
host_contigs="$(wc -l < "${HOST_NO_EBV}.fai")"
[[ "${host_contigs}" -eq "$((original_contigs - 1))" ]] || \
    die "PathSeq host reference did not remove exactly one contig"
grep -q "^${ebv_contig}"$'\t' "${HOST_NO_EBV}.fai" && die "${ebv_contig} remains in PathSeq host reference"

log "Downloading NCBI RefSeq viral genomes isolated from human hosts"
datasets download virus genome \
    taxon 10239 \
    --refseq \
    --host human \
    --include genome \
    --filename "${build_tmp}/viral_refseq.zip"
unzip -q "${build_tmp}/viral_refseq.zip" -d "${build_tmp}/viral_refseq"
viral_download="${build_tmp}/viral_refseq/ncbi_dataset/data/genomic.fna"
[[ -s "${viral_download}" ]] || die "NCBI data package did not contain genomic.fna"
mv "${viral_download}" "${MICROBE_FASTA}"
samtools faidx "${MICROBE_FASTA}"

grep -q '^NC_007605\.1'$'\t' "${MICROBE_FASTA}.fai" || \
    die "Viral reference does not contain EBV accession NC_007605.1"
grep -q '^NC_006273\.2'$'\t' "${MICROBE_FASTA}.fai" || \
    die "Viral reference does not contain CMV accession NC_006273.2"

log "Creating viral sequence dictionary"
gatk --java-options "-Xmx4g" CreateSequenceDictionary \
    --REFERENCE "${MICROBE_FASTA}" \
    --OUTPUT "${MICROBE_DICT}"

log "Downloading NCBI taxonomy and RefSeq catalog"
wget -q --retry-connrefused --waitretry=5 --tries=5 \
    -O "${TAXDUMP}" \
    "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
refseq_release="$(wget -qO- https://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER)"
[[ "${refseq_release}" =~ ^[0-9]+$ ]] || die "Invalid RefSeq release number: ${refseq_release}"
wget -q --retry-connrefused --waitretry=5 --tries=5 \
    -O "${REFSEQ_CATALOG}" \
    "https://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/RefSeq-release${refseq_release}.catalog.gz"

log "Building PathSeq taxonomy database"
gatk --java-options "-Xmx100g" PathSeqBuildReferenceTaxonomy \
    --reference "${MICROBE_FASTA}" \
    --output "${TAXONOMY_DB}" \
    --tax-dump "${TAXDUMP}" \
    --refseq-catalog "${REFSEQ_CATALOG}" \
    --min-non-virus-contig-length 0

log "Building viral BWA image"
gatk --java-options "-Xmx32g" BwaMemIndexImageCreator \
    --input "${MICROBE_FASTA}" \
    --output "${MICROBE_IMG}"

log "Building no-chrEBV host BWA image and k-mer hash"
gatk --java-options "-Xmx100g" BwaMemIndexImageCreator \
    --input "${HOST_NO_EBV}" \
    --output "${HOST_NO_EBV_IMG}"
gatk --java-options "-Xmx100g" PathSeqBuildKmers \
    --reference "${HOST_NO_EBV}" \
    --output "${HOST_NO_EBV_HSS}"

for output in \
    "${HOST_NO_EBV}" "${HOST_NO_EBV}.fai" "${HOST_NO_EBV_IMG}" "${HOST_NO_EBV_HSS}" \
    "${MICROBE_FASTA}" "${MICROBE_FASTA}.fai" "${MICROBE_DICT}" "${MICROBE_IMG}" \
    "${TAXDUMP}" "${REFSEQ_CATALOG}" "${TAXONOMY_DB}"; do
    [[ -s "${output}" ]] || die "Expected build output is missing or empty: ${output}"
done

rm -rf -- "${build_tmp}"
trap - EXIT