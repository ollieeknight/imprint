#!/usr/bin/env bash
# Run with bash or sbatch. The launcher submits one Slurm array task per donor.
#SBATCH --job-name=imprint-mgatk2-launch
#SBATCH --partition=compute
#SBATCH --exclude=s-sc-node038
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:15:00

set -euo pipefail

die() { echo "ERROR: $*" >&2; exit 1; }
abs() { python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "$1"; }

samplesheet=
cohort_dir=
fasta=
while (($#)); do
    case "$1" in
        --samplesheet)      samplesheet=${2:?}; shift 2 ;;
        --cohort-dir)       cohort_dir=${2:?}; shift 2 ;;
        --hardmasked-fasta) fasta=${2:?}; shift 2 ;;
        -h|--help)
            echo "Usage: bash $0 --samplesheet FILE --cohort-dir DIR --hardmasked-fasta FILE"
            exit 0
            ;;
        *) die "unknown argument: $1" ;;
    esac
done

[[ -f $samplesheet ]] || die "samplesheet not found: $samplesheet"
[[ -d $cohort_dir ]] || die "cohort directory not found: $cohort_dir"
[[ -f $fasta ]] || die "FASTA not found: $fasta"
[[ -f ${fasta}.fai ]] || die "FASTA index not found: ${fasta}.fai"
[[ -f ${fasta%.*}.dict ]] || die "FASTA dictionary not found: ${fasta%.*}.dict"
[[ -f ${fasta}.ann ]] || die "BWA-MEM3 index not found beside FASTA: ${fasta}.ann"

# Print the donor count, or the rows for one zero-based donor index.
read_sheet() {
    python3 - "$samplesheet" "$1" <<'PY'
import csv, os, re, sys

path, choice = sys.argv[1:]
required = ("patient", "cell_type", "status", "fastq_1", "fastq_2")
with open(path, newline="") as handle:
    reader = csv.DictReader(line for line in handle if line.strip())
    missing = [x for x in required if x not in (reader.fieldnames or [])]
    if missing:
        raise SystemExit("Missing samplesheet columns: " + ", ".join(missing))
    rows = [{k: (row.get(k) or "").strip() for k in required} for row in reader]

if not rows:
    raise SystemExit("Samplesheet contains no samples")
sample_status = {}
used_fastqs = set()
for line, row in enumerate(rows, 2):
    if not re.fullmatch(r"[A-Za-z0-9]+", row["patient"]):
        raise SystemExit(f"Row {line}: invalid patient")
    if not re.fullmatch(r"[A-Za-z0-9]+", row["cell_type"]):
        raise SystemExit(f"Row {line}: invalid cell_type")
    if row["status"] not in {"0", "1"}:
        raise SystemExit(f"Row {line}: status must be 0 or 1")
    if row["fastq_1"] == row["fastq_2"]:
        raise SystemExit(f"Row {line}: R1 and R2 are the same file")
    sample = row["patient"], row["cell_type"]
    if sample in sample_status and sample_status[sample] != row["status"]:
        raise SystemExit(f"Row {line}: sample has conflicting statuses")
    sample_status[sample] = row["status"]
    for key in ("fastq_1", "fastq_2"):
        if not os.path.isfile(row[key]):
            raise SystemExit(f"FASTQ not found: {row[key]}")
        if row[key] in used_fastqs:
            raise SystemExit(f"Row {line}: FASTQ is assigned more than once: {row[key]}")
        used_fastqs.add(row[key])

donors = sorted({row["patient"] for row in rows})
for donor in donors:
    normals = {row["cell_type"] for row in rows if row["patient"] == donor and row["status"] == "0"}
    if len(normals) != 1:
        raise SystemExit(f"Patient {donor}: expected one status=0 cell type, found {len(normals)}")

if choice == "count":
    print(len(donors))
else:
    donor = donors[int(choice)]
    for row in rows:
        if row["patient"] == donor:
            print("\t".join(row[x] for x in required))
PY
}

donor_count=$(read_sheet count)
samplesheet=$(abs "$samplesheet")
cohort_dir=$(abs "$cohort_dir")
fasta=$(abs "$fasta")

if [[ -z ${SLURM_ARRAY_TASK_ID:-} ]]; then
    mkdir -p "$cohort_dir/cohort/mgatk2_logs"
    exec sbatch \
        --job-name=imprint-mgatk2 \
        --partition=compute \
        --exclude=s-sc-node038 \
        --cpus-per-task=16 \
        --mem=64G \
        --time=48:00:00 \
        --array="0-$((donor_count - 1))" \
        --output="$cohort_dir/cohort/mgatk2_logs/%A_%a.log" \
        "$(abs "$0")" \
        --samplesheet "$samplesheet" \
        --cohort-dir "$cohort_dir" \
        --hardmasked-fasta "$fasta"
fi

CACHE=/sc-scratch/sc-scratch-cc12-ag-romagnani/apptainer_cache
ALIGN=$CACHE/align_process_env.sif
FASTP=$CACHE/quay.io-biocontainers-fastp-1.3.5--h43da1c4_0.img
SAMTOOLS=$CACHE/quay.io-biocontainers-samtools-1.23.1--ha83d96e_0.img
GATK=$CACHE/broadinstitute-gatk-4.6.2.0.img
MGATK2=$CACHE/quay.io-ollieeknight-mgatk2-2.0.0.img
BIND=(--bind /charite-store-f/f-cc12-ag-romagnani --bind /sc-projects/sc-proj-cc12-ag-romagnani --bind /sc-scratch/sc-scratch-cc12-ag-romagnani)

for image in "$ALIGN" "$FASTP" "$SAMTOOLS" "$GATK" "$MGATK2"; do
    [[ -f $image ]] || die "container not found: $image"
done

scratch=${SLURM_TMPDIR:-/sc-scratch/sc-scratch-cc12-ag-romagnani/nf_work_imprint}
work=$(mktemp -d "$scratch/imprint-mgatk2.XXXXXX")
trap 'rm -rf -- "$work"' EXIT
read_sheet "$SLURM_ARRAY_TASK_ID" > "$work/rows.tsv"
donor=$(awk 'NR == 1 {print $1}' "$work/rows.tsv")
mapfile -t cells < <(awk -F '\t' '!seen[$2]++ {print $2}' "$work/rows.tsv")

normal=
declare -A status bam
for cell in "${cells[@]}"; do
    status[$cell]=$(awk -F '\t' -v c="$cell" '$2 == c {print $3; exit}' "$work/rows.tsv")
    [[ ${status[$cell]} == 1 ]] || normal=$cell
done

# Realign bulk FASTQs to the NUMT-hardmasked reference using the steps in
# modules/align.nf. Keep the settings in sync if changing the nuclear pipeline:
# optical duplicate distance is fixed at 2500 here, with no front trimming.
# This script does not trim or tag duplex UMIs and must not be used for DupCaller runs.
echo "Processing $donor (normal: $normal)"
for cell in "${cells[@]}"; do
    sample=${donor}_${cell}
    sample_dir=$work/$sample
    mkdir -p "$sample_dir"
    lanes=()
    lane=0

    while IFS=$'\t' read -r _ row_cell _ r1 r2; do
        [[ $row_cell == "$cell" ]] || continue
        ((lane += 1))
        id=$(printf '%s_run%03d' "$sample" "$lane")
        prefix=$sample_dir/$id

        apptainer exec "${BIND[@]}" "$FASTP" fastp \
            --in1 "$r1" --in2 "$r2" \
            --out1 "${prefix}_R1_trimmed.fastq.gz" --out2 "${prefix}_R2_trimmed.fastq.gz" \
            --detect_adapter_for_pe --length_required 36 \
            --json "${prefix}_fastp.json" --html "${prefix}_fastp.html" --thread 8

        apptainer exec "${BIND[@]}" "$ALIGN" bwa-mem3 mem \
            -t 24 --bam=0 -Y -K 100000000 \
            -R "@RG\tID:${id}\tSM:${sample}\tPL:ILLUMINA\tLB:${sample}" \
            "$fasta" "${prefix}_R1_trimmed.fastq.gz" "${prefix}_R2_trimmed.fastq.gz" \
        | apptainer exec "${BIND[@]}" "$ALIGN" samtools fixmate \
            -@ 24 -m -O bam,level=1 - "${prefix}.fixmate.bam"

        apptainer exec "${BIND[@]}" "$SAMTOOLS" samtools sort \
            -@ 8 -m 5G -O bam -o "${prefix}.bam" "${prefix}.fixmate.bam"
        lanes+=("${prefix}.bam")
    done < "$work/rows.tsv"

    merged=$sample_dir/${sample}_merged_tagged.bam
    if ((${#lanes[@]} == 1)); then
        cp "${lanes[0]}" "$merged"
    else
        apptainer exec "${BIND[@]}" "$SAMTOOLS" samtools merge -f -@ 4 "$merged" "${lanes[@]}"
    fi
    apptainer exec "${BIND[@]}" "$SAMTOOLS" samtools index "$merged"

    dedup=$sample_dir/${sample}_dedup.bam
    apptainer exec "${BIND[@]}" "$GATK" gatk --java-options '-Xmx13107M' MarkDuplicates \
        --INPUT "$merged" --OUTPUT "$dedup" \
        --METRICS_FILE "$sample_dir/${sample}_markdup_metrics.txt" \
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT
    mv "$sample_dir/${sample}_dedup.bai" "${dedup}.bai"

    apptainer exec "${BIND[@]}" "$ALIGN" samtools sort -@ 2 -m 1G -n -O bam -o "$sample_dir/qname.bam" "$dedup"
    apptainer exec "${BIND[@]}" "$ALIGN" fgbio -Xmx4g --compression=0 --async-io=true CallOverlappingConsensusBases \
        --input "$sample_dir/qname.bam" --output "$sample_dir/clipped.qname.bam" \
        --metrics "$sample_dir/${sample}_overlap_metrics.txt" --ref "$fasta" --threads 2 \
        --agreement-strategy MaxQual --disagreement-strategy MaskLowerQual
    bam[$cell]=$sample_dir/${sample}_clipped.bam
    apptainer exec "${BIND[@]}" "$ALIGN" samtools sort -@ 2 -m 1G -O bam -o "${bam[$cell]}" "$sample_dir/clipped.qname.bam"
    apptainer exec "${BIND[@]}" "$ALIGN" samtools index "${bam[$cell]}"
done

for cell in "${cells[@]}"; do
    [[ ${status[$cell]} == 1 ]] || continue
    pair=${donor}_${cell}_v_${normal}
    result=$work/$pair
    mkdir "$result"
    apptainer exec "${BIND[@]}" "$MGATK2" mgatk2 paired \
        --tumor "${bam[$cell]}" --normal "${bam[$normal]}" \
        --reference "$fasta" --output "$result" --sample-name "$pair" --deduplication none

    # mgatk2 schema 3.0 emits the VCF, its index, and the callable BED. QC and
    # provenance are the ##mgatk2_qc= header line, not a sidecar JSON.
    for suffix in mt_variants.vcf.gz mt_variants.vcf.gz.tbi mt_callable.bed.gz; do
        [[ -s $result/$pair.$suffix ]] || die "missing MGATK2 output: $pair.$suffix"
    done

    destination=$cohort_dir/$donor/pairs/${cell}_v_${normal}/mgatk2
    incoming=${destination}.incoming.${SLURM_ARRAY_JOB_ID}
    mkdir -p "$(dirname "$destination")"
    rm -rf -- "$incoming"
    cp -a "$result" "$incoming"
    rm -rf -- "$destination"
    mv "$incoming" "$destination"
    echo "Installed $destination"
done
