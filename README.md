# imprint

imprint calls somatic variants in paired, FACS-sorted immune-cell populations from bulk WES, WGS or xGen UDSeq data. Each donor has one reference population and one or more populations to compare against it. The samplesheet calls these normal (`status=0`) and tumour (`status=1`), including in healthy donors. The normal population may share clonal variants with the tumour population; it is not necessarily a constitutional germline control.

This pipeline uses Nextflow DSL2 and was initially designed to run on the Charité SC HPC. Paths and resource settings in this repository are for that cluster and need changing for another site.

## Requirements

- Nextflow 24.10 or newer, SLURM and Apptainer. Bulk calling also needs a GPU for DeepSomatic.
- Paired-end FASTQs and a samplesheet as described below.
- GRCh38 references, indexes, annotation resources and container images configured in [conf/resources.config](conf/resources.config). Capture data also need matching intervals in [conf/probekits.config](conf/probekits.config).

Set the work and cache directories in `conf/resources.config`, and the queues, bind paths and resource limits in [conf/process.config](conf/process.config). Build the shared BWA-MEM3 alignment image from [assets/containers/align_env.def](assets/containers/align_env.def). Reference preparation scripts are in [assets/reference](assets/reference).

DupCaller needs its own image and reference indexes; see below. Optional analyses need only the resources for the extras you select.

## Samplesheet

Use a CSV with these columns:

```csv
patient,cell_type,status,fastq_1,fastq_2
HC01,CD14,0,/data/HC01_CD14_R1.fastq.gz,/data/HC01_CD14_R2.fastq.gz
HC01,NKG2A,1,/data/HC01_NKG2A_R1.fastq.gz,/data/HC01_NKG2A_R2.fastq.gz
HC01,NKG2C,1,/data/HC01_NKG2C_R1.fastq.gz,/data/HC01_NKG2C_R2.fastq.gz
```

`patient` identifies the donor; `cell_type` identifies the population. Both must contain only letters and numbers. Each donor needs exactly one distinct normal sample. Repeat rows for additional sequencing lanes of the same sample; imprint merges them before duplicate marking. Use absolute FASTQ paths, with each file appearing only once.

## Running

For bulk WES:

```bash
nextflow run /path/to/imprint/main.nf \
  -profile slurm,xgen_exome_v2 \
  --samplesheet /path/to/samplesheet.csv \
  --cohort_name my_cohort \
  --outdir /path/to/outdir
```

Choose the profile that matches the capture kit: `agilent_v6`, `agilent_v7`, `agilent_v8`, `twist_v2` or `xgen_exome_v2`. Use `-profile slurm,wgs` for WGS, or `-profile slurm,xgen_exome_v2,dupcaller` for xGen UDSeq. Do not combine WGS and capture profiles. Add `-resume` to reuse completed tasks after fixing a failed run.

Bulk calling uses Mutect2, Strelka2 and DeepSomatic. Manta supplies candidate indels to Strelka2 and reports structural variants. The small-variant union retains sites that PASS any caller, then adds Vafator read evidence, VEP annotations and VarLociraptor probabilities. Raw and filtered caller outputs remain available. Population-frequency annotations do not impose a final study filter.

Defaults are in [nextflow.config](nextflow.config). WES uses 50 bp padded capture intervals by default; `--off_target false` restricts calling to the targets. BQSR is off for WES and on for WGS. Set `--optical_dup_dist` for the sequencing instrument: the default is 2500 for patterned flow cells; unpatterned instruments use 100. `--trim_front` removes a fixed number of bases from both reads in bulk mode only.

## xGen UDSeq

The DupCaller profile accepts xGen Exome Hyb Panel v2 libraries with an 8 bp UMI at the start of each mate. imprint corrects each UMI against the [xGen UMI32 allowlist](assets/umi_allowlists/xgen_8bp_umi32.txt), allowing one mismatch, then records the corrected pair in the read name and `DB` tag. It checks the tags after duplex-aware duplicate marking and after the CRAM round trip.

DupCaller uses raw read families, without molecular consensus reads, overlap correction or BQSR. Each pair runs `call` and `estimate`, followed by a cohort `summarize`. VEP annotates PASS SBS and indel calls while preserving the DupCaller family fields. The bulk callers, Vafator, VarLociraptor and Mutect2 panels of normals are excluded from this mode.

The image definition pins DupCaller dev commit `53eb785` (version 1.1.2):

```bash
apptainer build dupcaller_1.1.2-dev53eb785.sif assets/containers/dupcaller_env.def
```

Point `container_dupcaller` to that image. Keep the version and commit in its filename; imprint reads the version for provenance. DupCaller reuses the bulk FASTA and BWA-MEM3 index and needs five additional files beside the FASTA: `<fasta>.ref.h5`, `<fasta>.tn.h5`, `<fasta>.hp.h5`, `<fasta>.str.h5` and `<fasta>.dbs.h5`.

[bin/build_dupcaller_refs.sh](bin/build_dupcaller_refs.sh) prepares these indexes. The repeat scan needs PERF installed separately. It uses `PERF.core -m 1 -M 10 -u 2`; keep `-u 2` so short, low-copy repeats are included. The index step passes that repeat TSV to `DupCaller.py index -rt`.

Configure noise masks and an optional indel ePoN in `conf/resources.config`, with a matching `.tbi` for each file. Check that masks match the chemistry, reference build and contig names. Unless overridden, the maximum zero-quality fraction is 0.5 with masks and 0.1 without them. Set `--dupcaller_seed` to reproduce the detection-power simulation; otherwise DupCaller records its chosen seed in the call log.

## Outputs

Within the output directory:

```text
cohort/
  samplesheet.csv
  run_params.json
  output_manifest.json
  multiqc/
  somalier/
  dupcaller/                         # DupCaller cohort summary
{donor}/
  samples/{cell_type}/
    alignment/                      # CRAM and index
    qc/
  pairs/{tumour}_v_{normal}/
    variant_calling/                 # bulk calls, including raw/
    dupcaller/                      # DupCaller runs
      calls/
      burden/
      annotated/                    # VEP-annotated PASS SBS and indels
```

`output_manifest.json` (schema 2.0) lists the emitted files and their sample or pair metadata. `run_params.json` records effective parameters, references, containers and launch details. Bulk final calls are `{pair_id}.somatic.vcf.gz` under `variant_calling/`.

DupCaller keeps its `SBS/`, `INDEL/`, `DBS/` and `ERROR/` subdirectories. Callsets include the original VCF and a sorted, bgzip-compressed copy with a tabix index. The uncompressed `_fail.vcf` files retain rejected candidates and their filter reasons. `burden/` contains the estimates and plots produced for that callset.

QC includes fastp, Mosdepth, VerifyBamID2, Somalier, Riker and MultiQC. In DupCaller runs, distinguish Mosdepth read depth from molecular duplex coverage in the DupCaller coverage BED. DupCaller error profiles replace the bulk read-position error report.

## Optional analyses

Add a comma-separated list, for example `--extras hla,kir,telseq`:

| Extra | Analysis | Additional requirement |
| --- | --- | --- |
| `hla` | OptiType HLA typing per donor | HLA reference and Yara index |
| `kir` | kir-mapper copy number and genotyping per donor | kir-mapper image and database |
| `pathseq` | PathSeq microbial screening per sample | Host, microbial and taxonomy references |
| `mixcr` | MiXCR immune-receptor repertoire | MiXCR licence (`--mixcr_license`) |
| `telseq` | TelSeq telomere-length estimate | TelSeq image |

These analyses may fail without cancelling the core outputs. MiXCR uses the trimmed reads for the selected mode; TelSeq uses merged BAMs before duplicate marking.

Mitochondrial calling runs separately through [bin/imprint-mgatk2.sbatch](bin/imprint-mgatk2.sbatch), for bulk data only. It needs a NUMT-hardmasked reference and mgatk2 2.0.0 or newer. The launcher submits one SLURM array task per donor. Each pair produces `.mt_variants.vcf.gz`, its `.tbi` and `.mt_callable.bed.gz` under `{donor}/pairs/{tumour}_v_{normal}/mgatk2/`. VCF samples are ordered `NORMAL`, `TUMOR`; the `##mgatk2_qc=` header contains QC and provenance. Use the callable BED as the denominator for mitochondrial burden. These outputs are outside the Nextflow manifest and must be found by path.

## Local checks

```bash
python3 bin/trim_dupcaller_umis.py --self-test
python3 -m py_compile bin/trim_dupcaller_umis.py bin/collate_kir.py assets/reference/prepare_probe_kits.py
nextflow lint main.nf modules/*.nf subworkflows/*.nf
```

These checks do not validate the cluster containers, reference files or biological results. Test those on the cluster before using a release for analysis.

## Citation and licence

imprint is released under the MIT licence; see [LICENSE](LICENSE). If you use it
in published work, cite it using the metadata in [CITATION.cff](CITATION.cff).
