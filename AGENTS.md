# imprint agent guide

Read this guide before editing and keep user instructions in `README.md`.

## Safety and scope

- Production runs on the Charité SC HPC, using SLURM with Apptainer.
- Never commit, push, or do remote writes.
- Keep mitochondrial calling outside the nuclear Nextflow graph, with bulk input only.

## Biological model

Samplesheet `status=1` means tumour; `status=0` means matched normal. The normal is not necessarily a constitutional germline control and may share biological signal with the tumour.

Two primary paths:

```text
bulk WES/WGS
  -> fastp -> BWA-MEM3 -> merge -> MarkDuplicates -> overlap correction
  -> Mutect2 + Manta/Strelka2 + DeepSomatic
  -> PASS union -> Vafator/VEP + VarLociraptor

xGen UDSeq
  -> allowlist-corrected 8+8 bp UMI trim -> BWA-MEM3 -C
  -> merge -> GATK duplex-aware MarkDuplicates -> tag validation
  -> DupCaller call -> estimate -> cohort summarize
```

DupCaller uses raw read families without molecular consensus reads. Keep its outputs out of the bulk callers, Manta, BQSR, Vafator, VarLociraptor and Mutect2 PoN. Both modes share `VEP_ANNOTATE`: DupCaller PASS SBS and indel VCFs are published to `dupcaller/annotated/`; bulk VEP output is an unpublished intermediate for `VARLOCIRAPTOR_MERGE`. The process takes `mutation_type` to name its output and uses `publishDir ..., enabled: params.dupcaller` to control publication.

imprint tracks the DupCaller dev branch. Build the image from `assets/containers/dupcaller_env.def`, pinned to a commit. Include the version and commit in `container_dupcaller`. `main.nf` reads the version from the tag, which must agree with `src/DupCaller_sub/__init__.py`. Revisit the pin when dev merges into main.

## Hard invariants

- Samplesheet columns: `patient,cell_type,status,fastq_1,fastq_2`.
- `patient` and `cell_type` must be non-empty `[A-Za-z0-9]+`.
- `status` = `0` or `1`; every donor has exactly one distinct reference sample.
- Repeated rows for one donor/cell type/status = lanes; duplicate FASTQ paths or pairs fail.
- WGS and capture profiles cannot mix.
- DupCaller requires `xgen_exome_v2`; reject WGS, Agilent, Twist and mixed chemistry before scheduling.
- Explicit resources must exist; indexed resources must have matching `.tbi` files.
- DupCaller needs five h5 files beside the FASTA: `ref`, `tn`, `hp`, `str`, `dbs`. All five use the FASTA filename as their prefix. Build `str` from a PERF repeat TSV passed to `DupCaller.py index -rt`.
- Reject unknown extras.

## DupCaller chemistry and BAM contract

Authoritative allowlist: `assets/umi_allowlists/xgen_8bp_umi32.txt`.

- R1 and R2 each begin with one 8 bp UMI.
- Correct each component independently by at most one mismatch.
- Query-name suffix = `_R1UMI+R2UMI`; FASTQ/BAM tag = `DB:Z:R1UMI-R2UMI`.
- Query-name and `DB` values must stay equal after alignment, merge, duplicate marking, CRAM round trip.
- Invalid, short, unmatched, unequal, malformed, or mate-mismatched FASTQ input cannot silently enter a family.
- Preserve RG, `AS`, `XS`, `NM`, mate fields, CIGAR, flags, duplicate flags, GATK `DT` behaviour.
- Do not correct overlaps, cap qualities, recalibrate bases or realign reads before DupCaller.

`bin/trim_dupcaller_umis.py` handles barcode correction. `DUPCALLER_VALIDATE_BAM` enforces the post-alignment contract. Keep one implementation of each rule.

On dev, `--rescue` is `store_true`; passing a value aborts the run. `-gaf` takes a float, supplied through `dupcaller_germline_af_cutoff`.

## Bulk evidence rules

Keep bulk behaviour independent:

- Mutect2 must use zero candidate LOD settings; never add a project VAF floor.
- Keep `FilterMutectCalls` and both raw/filtered VCFs.
- Publish Manta SVs and pass its candidate indels to Strelka2.
- Retain PASS calls from any bulk caller in the all-PASS union.
- Use VEP gnomAD tracks for annotation only.
- `TUMOR_NORMAL_FILTER` marks territory and passes records through; downstream study analyses apply filters.
- Keep the external Mutect2 PoN and run-derived cohort PoN separate and bulk-only.

## Ownership

| Path | Owner |
| --- | --- |
| `main.nf` | startup checks, metadata, top graph, provenance, manifest records |
| `nextflow.config` | general user parameters |
| `conf/probekits.config` | assay profiles and capture intervals |
| `conf/resources.config` | HPC paths, references, containers |
| `conf/process.config` | executors, retries, resources |
| `modules/*.nf` | process command, container, inputs, outputs, publish policy |
| `subworkflows/*.nf` | channel orchestration |
| `modules/reporting.nf` | manifest/provenance schema |
| `bin/trim_dupcaller_umis.py` | xGen UMI correction and strict FASTQ pairing |

Trace output changes through the producer, subworkflow, `main.nf`, manifest, downstream readers and README. Manifest records come from emitted channels, never predicted files.

`bin/imprint-mgatk2.sbatch` runs under SLURM, outside Nextflow. It emits no
channel, so its outputs are absent from `output_manifest.json` even though
they sit inside `{donor}/pairs/`. Readers must find mitochondrial calls by
path. Do not add predicted mgatk2 paths to the manifest. Manifest coverage
would require integrating mgatk2 into the Nextflow graph.

## Metadata contract

Sample metadata, built once in `main.nf` from samplesheet:

| Key | Meaning |
| --- | --- |
| `id` | `{donor}_{cell_type}`, and read-group `SM` |
| `donor` | samplesheet `patient` |
| `cell_type` | samplesheet `cell_type`, and sample output directory |
| `status` | `0` normal, `1` tumour |
| `donor_sample_count` | donor cardinality for `groupKey` |
| `run_count`, `run_id` | lane fields; `sampleMeta()` removes them at merge to keep them out of downstream task hashes. Retain lane identity only in lane filenames and read-group IDs. |
| `short_inserts`, `read_length` | fastp-derived, added to tumour samples before pairing |

Pair metadata, built once where tumour and normal meet:

| Key | Meaning |
| --- | --- |
| `pair_id` | `{donor}_{pair_dir}` |
| `pair_dir` | `{tumor_cell_type}_v_{normal_cell_type}`, pair output directory |
| `donor` | donor identifier |
| `tumor_id`, `normal_id` | sample IDs; names also match manifest schema 2.0 |
| `tumor_cell_type`, `normal_cell_type` | cell types; processes must read these fields rather than parse identifiers |
| `short_inserts` | tumour soft-clip autodetection for Mutect2 |

Read cell types and output directories from metadata; do not reconstruct
them by stripping identifier prefixes. Add a metadata field where needed.

Use the tag for the process scope:

| Scope | Tag |
| --- | --- |
| cohort | `'cohort'`, or `'intervals'` for the two interval builders |
| donor | `"${meta.donor}"` |
| sample | `"${meta.id}"`, or `"${sample_id}"` for a per-sample stream inside a pair |
| pair | `"${meta.pair_id}"` |

## Channel rules

- Each biological tuple starts with metadata: `[meta, ...]`.
- Preserve metadata during grouping.
- Per-sample streams inside a pair carry `[meta, role, sample_id, ...]`, where `role` = `tumour` or `normal`. Branch on `role`, never on the status flag.
- Join biological data by sample ID, pair ID or donor, never by incidental filenames.
- Use `groupKey` with known `run_count`, interval count, or donor sample count.
- Use strict joins for required complete sets.
- Keep one base resource label per process; modifiers may add runtime or memory behaviour.
- Name subworkflow emits unless there is only one; `nextflow lint` rejects a lone named emit.

## Output contracts

Common:

```text
cohort/{samplesheet.csv,run_params.json,output_manifest.json}
cohort/{multiqc,somalier}/
{donor}/samples/{cell_type}/{alignment,qc}/
```

Mitochondrial pair outputs live under `{donor}/pairs/{tumor}_v_{normal}/mgatk2/`
and are exactly three files: `.mt_variants.vcf.gz`, its `.tbi`, and
`.mt_callable.bed.gz`. Do not publish a TSV, sidecar JSON or log in the pair
directory. mgatk2 schema 3.0 stores QC/provenance in the `##mgatk2_qc=` VCF header. Sample order in that VCF is
`NORMAL`, `TUMOR`. The suffix check in the sbatch is the contract; update both
together.

Bulk pair calls stay under `variant_calling/`. DupCaller pair outputs live under:

```text
{donor}/pairs/{tumor}_v_{normal}/dupcaller/calls/
{donor}/pairs/{tumor}_v_{normal}/dupcaller/burden/
cohort/dupcaller/
```

Publish `calls/` and `burden/` as whole directories, preserving DupCaller's `SBS/`, `INDEL/`, `DBS/` and `ERROR/` subfolders. Build manifest records by walking the emitted directory. Derive `kind` from the filename alone, which already includes the mutation type. Preserve the directory layout: `estimate` and `summarize` read it back.

Manifest schema 2.0 records `library_mode`. DupCaller runs also record the version and `xgen_udseq_8bp_umi32`; derive the version from `container_dupcaller`. Preserve raw DupCaller VCF family INFO/FORMAT fields in downstream transformations.

## Optional analyses

Accepted `--extras` are `hla,kir,pathseq,mixcr,telseq`, with existing aliases. These provide characterisation, not primary-call evidence. Each must use `errorStrategy ignore` on its final attempt so it cannot cancel core outputs. Preserve that policy for new extras. MiXCR consumes mode-appropriate trimmed FASTQs; TelSeq consumes merged BAMs before duplicate marking. Keep `bin/imprint-mgatk2.sbatch` bulk-only. mgatk2 >= 2.0.0 (schema 3.0) uses `--tumor`/`--normal` and produces only the VCF, `.tbi` and callable BED, with QC in the `##mgatk2_qc=` header.

## Nextflow traps

- Staged inputs are usually symlinks; write new outputs rather than edit staged paths.
- Exit 127 = command/container mismatch; split container boundaries rather than inflate unrelated images.
- Container `site-packages` is read-only at runtime. Packages that cache in their installation directory fail with `Errno 30`; redirect the cache to the task directory. For `sigProfilerPlotting`, use `SIGPROFILERPLOTTING_VOLUME`.
- Check the process exit status before investigating shutdown messages.
- Corrected process scripts rerun under `-resume`; valid upstream work stays cached.
- Output declarations must match actual pinned DupCaller image files, including empty callsets. `DUPCALLER_ESTIMATE` identifies new files rather than hard-coding a list because burden outputs vary with the callset and SigProfilerPlotting naming.
- `publishDir` must use closure form (`publishDir { "..." }`) to see task inputs. A bare string resolves in script scope and fails with `No such variable: meta`. `tag` tolerates a bare string; `publishDir` does not.
- A process that declares both a directory output and a glob inside it publishes neither. Declare the directory alone and exclude unwanted files by returning `null` from `saveAs`.
- A trailing `[ -f x ] && mv x y` sets the script exit status; a missing file fails the task.
- GATK `-L` accepts `.bed`, `.interval_list`, `.intervals` and VCF, but not a sequence dictionary.
- `file(null)` throws before any custom message, so null-check path parameter before `file()`.

## Completion checks

1. Inspect `git status --short`; account for every touched file.
2. Run `python3 bin/trim_dupcaller_umis.py --self-test` and Python compile when Python changes.
3. Run `bash -n` on changed shell scripts.
4. Parse every R file without executing analyses.
5. Run `nextflow lint main.nf modules/*.nf subworkflows/*.nf`.
6. Search for removed parameters/processes, stale output paths, undeclared parameters, orphan manifest kinds.
7. Run `git diff --check` and review complete diff.
8. State that HPC containers/references and biological execution remain unverified until tested on the cluster.
