#!/usr/bin/env bash
#SBATCH --job-name=imprint-mito-ref
#SBATCH --partition=compute
#SBATCH --exclude=s-sc-node038
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=/sc-projects/sc-proj-cc12-ag-romagnani/ref/imprint/genome/mito-hardmask.%j.log
#SBATCH --error=/sc-projects/sc-proj-cc12-ag-romagnani/ref/imprint/genome/mito-hardmask.%j.log

set -euo pipefail

REF=/sc-projects/sc-proj-cc12-ag-romagnani/ref/imprint/genome/GRCh38_no_alt_analysis_set.fasta
OUTDIR=/sc-projects/sc-proj-cc12-ag-romagnani/ref/imprint/genome/mito-hardmasked
FASTA=$OUTDIR/GRCh38_no_alt_analysis_set.hardmasked.fasta
DICT=$OUTDIR/GRCh38_no_alt_analysis_set.hardmasked.dict
BLACKLIST=$OUTDIR/hg38.full.blacklist.bed

BEDTOOLS=/sc-scratch/sc-scratch-cc12-ag-romagnani/apptainer_cache/quay.io-biocontainers-bedtools_2.27.1--h077b44d_9.img
SAMTOOLS=/sc-scratch/sc-scratch-cc12-ag-romagnani/apptainer_cache/quay.io-biocontainers-samtools-1.23.1--ha83d96e_0.img
ALIGN=/sc-scratch/sc-scratch-cc12-ag-romagnani/apptainer_cache/align_process_env.sif
BIND=(--bind /sc-projects/sc-proj-cc12-ag-romagnani --bind /sc-scratch/sc-scratch-cc12-ag-romagnani)

mkdir -p "$OUTDIR"
rm -f "$FASTA" "$FASTA".* "$DICT"

curl -fsSL \
    https://raw.githubusercontent.com/caleblareau/mitoblacklist/master/combinedBlacklist/hg38.full.blacklist.bed \
    -o "$BLACKLIST"

apptainer exec "${BIND[@]}" "$BEDTOOLS" \
    bedtools maskfasta -fi "$REF" -bed "$BLACKLIST" -fo "$FASTA"

apptainer exec "${BIND[@]}" "$SAMTOOLS" samtools faidx "$FASTA"
apptainer exec "${BIND[@]}" "$SAMTOOLS" samtools dict -o "$DICT" "$FASTA"

THREADS="${SLURM_CPUS_PER_TASK:-8}"

apptainer exec "${BIND[@]}" "$ALIGN" bwa-mem3 index -t "$THREADS" --max-memory 80G "$FASTA"
