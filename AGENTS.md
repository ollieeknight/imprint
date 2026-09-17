# AGENTS.md

imprint is a Nextflow DSL2 somatic-calling pipeline. Read [README.md](README.md) for what it does, how to run it and what it emits; this file covers only what the code does not confess.

## Before you claim a change works

```bash
nextflow lint main.nf modules/*.nf subworkflows/*.nf
python3 -m py_compile bin/trim_dupcaller_umis.py bin/collate_kir.py assets/reference/prepare_probe_kits.py
python3 bin/trim_dupcaller_umis.py --self-test
```

That is the whole local gate. It parses code; it proves nothing about containers, references or biological output. Real runs happen on the Charité SC HPC under SLURM with Apptainer, and you cannot reach it from here, so report a change as "lints clean, untested on cluster" and let a human run it.

## Layout

`main.nf` holds samplesheet preflight, channel construction and pair formation. `subworkflows/*.nf` wire processes into a stage; `modules/*.nf` define the processes themselves. Each subworkflow has a same-named module file. Config splits three ways: `nextflow.config` for params and defaults, `conf/resources.config` for site paths and container images, `conf/process.config` for queues, binds and label-to-resource mapping, `conf/probekits.config` for capture-kit profiles.

## Invariants

- **`meta` is the contract.** Sample-level channels carry `id`, `donor`, `cell_type`, `status`; pair-level channels carry `pair_id`, `pair_dir`, `tumor_id`, `normal_id`. Adding a key to a meta map in one process without threading it through the joins upstream silently drops records in `join`/`combine`.
- **Publish paths are the output schema.** Every `publishDir` closure follows `${params.outdir}/${meta.donor}/samples/${meta.cell_type}/…` or `…/pairs/${meta.pair_dir}/…`. Changing one changes the documented output tree, `output_manifest.json` (schema 2.0) and downstream consumers.
- **Resources come from labels, not from numbers in the module.** Use the existing `process_*` labels; add a per-process `memory`/`cpus` closure only when it scales off input size, as the alignment processes do.
- **Filenames are derived, not passed.** Processes compute output names by `.replace()` on the input filename so `-resume` stays stable. Keep that idiom rather than introducing separate name parameters.
- **Two calling modes share one pipeline.** Bulk (Mutect2/Strelka2/DeepSomatic + Vafator + VarLociraptor) and DupCaller raw-family calling are mutually exclusive; DupCaller excludes the bulk callers, PoN and VarLociraptor entirely. A change to shared alignment or QC must hold for both.
- **Normal is not germline.** `status=0` is the donor's reference population and can share clonal variants with the tumour population. Never introduce logic that assumes the normal is a constitutional control.
- **Optional extras fail soft.** An `--extras` analysis that dies must not cancel core outputs.

## Conventions

Container images are always `params.container_*`, never hardcoded. Version pins that carry provenance (DupCaller's `1.1.2-dev53eb785`) live in image filenames and are parsed at runtime — keep the format. Params get a comment only when the value is non-obvious (why 2500, why `-u 2`); the name carries the rest.

## Git

`main` is the default branch; never git add, commit or push without human review.
