# imprint

imprint calls somatic variants between FACS-sorted immune-cell populations from the same donor, using bulk WES, WGS or xGen UDSeq duplex data. Each donor has one reference population (`status=0`, called the normal) and one or more populations compared against it (`status=1`, the tumour), healthy donors included. The normal is a cell population, not a germline control, and can share clonal variants with the tumour.

It is a Nextflow pipeline built for the Charité SC HPC (SLURM + Apptainer). Paths and resources in `conf/` are for that cluster.

## Quick start

Write a samplesheet, one row per lane:

```csv
patient,cell_type,status,fastq_1,fastq_2
HC01,CD14,0,/data/HC01_CD14_R1.fastq.gz,/data/HC01_CD14_R2.fastq.gz
HC01,NKG2A,1,/data/HC01_NKG2A_R1.fastq.gz,/data/HC01_NKG2A_R2.fastq.gz
HC01,NKG2C,1,/data/HC01_NKG2C_R1.fastq.gz,/data/HC01_NKG2C_R2.fastq.gz
```

`patient` and `cell_type` must be alphanumeric, each donor needs exactly one normal, and each FASTQ may appear only once. Rows that share `patient` and `cell_type` are lanes of one sample and are merged.

Then run:

```bash
nextflow run /path/to/imprint/main.nf \
  -profile slurm,xgen_exome_v2 \
  --samplesheet samplesheet.csv \
  --cohort_name my_cohort \
  --outdir results
```

Add `-resume` to pick up a failed run where it stopped.

## Modes

| Profile | Data | Calling |
| --- | --- | --- |
| `slurm,<kit>` | WES (`agilent_v6`, `agilent_v7`, `agilent_v8`, `twist_v2`, `xgen_exome_v2`) | Mutect2, Strelka2 and DeepSomatic, merged |
| `slurm,wgs` | WGS | Same, genome-wide, with BQSR |
| `slurm,xgen_exome_v2,dupcaller` | xGen UDSeq | DupCaller on raw duplex families |

**Bulk** keeps every site that any caller passes, then adds Vafator read counts, VEP annotation and VarLociraptor probabilities. No population-frequency filter is applied; that is left to the analysis. Exomes are called on 50 bp padded targets unless you pass `--off_target false`. Bulk calling needs a GPU for DeepSomatic.

**DupCaller** trims and error-corrects the 8 bp UMIs against the [UMI32 allowlist](assets/umi_allowlists/xgen_8bp_umi32.txt), marks duplicates duplex-aware, and runs DupCaller `call`, `estimate` and a cohort `summarize`. No consensus reads, overlap correction or BQSR, since DupCaller models raw families itself. VEP annotates the SBS and indel VCFs; select `FILTER=PASS` for final calls.

Set `--optical_dup_dist 100` for unpatterned flow cells (HiSeq 2500, MiSeq, NextSeq 550). Other defaults are in [nextflow.config](nextflow.config).

## Setup

All site settings live in `conf/`: references and container images in [resources.config](conf/resources.config), queues, binds and resources in [process.config](conf/process.config), capture kits in [probekits.config](conf/probekits.config).

Build the two local images:

```bash
apptainer build align_process_env.sif assets/containers/align_env.def
apptainer build dupcaller.sif assets/containers/dupcaller_env.def   # UDSeq only
```

Reference preparation scripts are in [assets/reference](assets/reference). DupCaller also needs five index files beside the FASTA (`<fasta>.ref.h5`, `.tn.h5`, `.hp.h5`, `.str.h5`, `.dbs.h5`); [bin/build_dupcaller_refs.sh](bin/build_dupcaller_refs.sh) builds them and needs PERF installed. Noise masks and an indel ePoN are set in `resources.config` and should match the chemistry and reference build.

## Outputs

```text
cohort/
  run_params.json            # parameters, references, containers
  output_manifest.json       # every emitted file with its sample/pair metadata
  multiqc/  somalier/  dupcaller/
{donor}/
  samples/{cell_type}/
    alignment/               # CRAM
    qc/
  pairs/{tumour}_v_{normal}/
    variant_calling/         # {pair_id}.somatic.vcf.gz; per-caller files in raw/
    dupcaller/               # calls/, burden/, annotated/
```

QC is fastp, Mosdepth, VerifyBamID2, Somalier and Riker, collected in MultiQC. In DupCaller runs, use the coverage BED in `dupcaller/calls/` for duplex depth; Mosdepth reports deduplicated read depth.

## Optional analyses

Add any of these with `--extras`, e.g. `--extras hla,kir,telseq`. A failing extra does not stop the core pipeline.

| Extra | Runs | Needs |
| --- | --- | --- |
| `hla` | OptiType HLA typing, per donor | HLA reference with its Yara index |
| `kir` | kir-mapper copy number and genotypes, per donor | kir-mapper image and database |
| `pathseq` | PathSeq microbial screen, per sample | Host, microbe and taxonomy references |
| `mixcr` | MiXCR immune-receptor repertoire | MiXCR licence (`--mixcr_license`) |
| `telseq` | TelSeq telomere length | TelSeq image |

## Mitochondrial variants

These run outside Nextflow, for bulk data only: [bin/imprint-mgatk2.sbatch](bin/imprint-mgatk2.sbatch) submits one SLURM job per donor and writes `{donor}/pairs/{pair}/mgatk2/`. It needs a NUMT-hardmasked reference and mgatk2 2.0.0 or newer. Use `*.mt_callable.bed.gz` as the burden denominator. These files are not in `output_manifest.json`.

## Citation

imprint is MIT-licensed. If you use it in published work, cite it with [CITATION.cff](CITATION.cff).
