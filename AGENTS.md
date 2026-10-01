# AGENTS.md

imprint is a Nextflow DSL2 somatic-calling pipeline. [README.md](README.md) covers what it does, how to run it and where config lives; this file covers what the code doesn't show.

## Before you claim a change works

```bash
nextflow lint main.nf modules/*.nf subworkflows/*.nf
python3 -m py_compile bin/trim_dupcaller_umis.py bin/collate_kir.py assets/reference/prepare_probe_kits.py
python3 bin/trim_dupcaller_umis.py --self-test
bash bin/build_dupcaller_refs.sh --self-test
```

That is the whole local gate. It parses code and proves nothing about containers, references or biological output. Real runs happen on the Charité SC HPC, which you cannot reach, so report a change as "lints clean, untested on cluster" and let a human run it.

## Layout

`main.nf` does samplesheet preflight, channel construction and pair formation. `modules/*.nf` define processes; `subworkflows/*.nf` wire them into stages, include from any module and can nest (`somatic_calling` wraps the per-caller subworkflows), so follow the `include` lines rather than assuming file names pair up. Shared helpers (reference inputs, the calling BED, once-per-normal grouping) live in `subworkflows/common.nf`.

## Invariants

- **`meta` is the contract.** Sample channels carry `id`, `donor`, `cell_type`, `status`; pair channels carry `pair_id`, `pair_dir`, `tumor_id`, `normal_id`. A key added in one place but not threaded through upstream joins silently drops records in `join`/`combine`.
- **Publish paths are the output schema.** Every `publishDir` follows `${params.outdir}/${meta.donor}/samples/${meta.cell_type}/…` or `…/pairs/${meta.pair_dir}/…`. Changing one changes the documented output tree, `output_manifest.json` and everything downstream.
- **`-resume` must stay stable.** Output names are derived by `.replace()` on input names, never passed in. Anything gathered with `groupTuple` or `collect` arrives in task-completion order, so sort it by file name before it reaches a process; otherwise the command changes and every downstream task reruns.
- **Resources come from labels.** Use the existing `process_*` labels; add a per-process `memory`/`cpus` closure only when it scales off input size, as alignment does.
- **Two calling modes share one pipeline.** Bulk (Mutect2/Strelka2/DeepSomatic, Vafator, VarLociraptor, PoN) and DupCaller are mutually exclusive. Changes to shared alignment or QC must hold for both.
- **Normal is not germline.** `status=0` is the donor's reference population and can share clonal variants with the tumor, so no logic may treat it as a constitutional control.
- **Optional extras fail soft.** A failing `--extras` analysis must not cancel core outputs.

## Conventions

Container images are always `params.container_*`, never hardcoded. A channel can feed several consumers directly; don't `multiMap` just to copy it. Comment a param only when its value is non-obvious (why `optical_dup_dist = 2500`); the name carries the rest.

## Git

`main` is the default branch. Never git add, commit or push without human review.
