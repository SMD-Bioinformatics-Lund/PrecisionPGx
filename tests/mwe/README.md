# MWE harnesses for issues #26 and #27

Two tiny Nextflow scripts and a runner that demonstrate the bugs *and*
the fixes side-by-side. Used to verify the changes on the
`fix-pharmcat-pairing-and-staging` branch.

## What the bugs are

### #26 — PharmCAT re-indexes the reference genome every time

The pipeline passes `reference.fna.bgz` and its `.fai` into
`PHARMCAT_VCFPREPROCESSOR`. With Nextflow's default `symlink` staging,
the workdir's `reference.fna.bgz.fai` is a *symlink to the shared source*.
htslib, on opening the bgzipped fasta, decides the `.fai` is stale (or
just touches it on close) and writes — and that write reaches the
*source* file through the symlink. Under any parallel run, multiple
tasks race on the source-side `.fai`, producing the `Could not
understand FASTA index` errors reported in the issue.

**Fix:** [`stageInMode = 'copy'`](../../conf/modules/pharmcat_vcf_processing.config) on
`PHARMCAT_VCFPREPROCESSOR`. Inputs are copied into each task's workdir,
so htslib's write hits a private copy. The shared source is untouched
regardless of how many tasks run in parallel.

### #27 — BCFTOOLS_VIEW filters with the wrong sample's BED

The VCF stream and the target-pass BED stream feeding BCFTOOLS_VIEW
come from different upstream subworkflows (PHARMCAT_VCFPREPROCESSOR
and QC_BAM). Their meta dicts share `id` but diverge in other fields
(e.g. `data_type`). The original code stripped meta from the BED
channel — `.map { meta, bed -> [bed] }` — so Nextflow paired the two
channels by emission order. Under genuinely parallel execution the
orderings diverge and the wrong BED reaches the wrong sample, silently.

A naive switch to `.join(other, ...)` doesn't help either: comparing on
the full meta map fails to match because of the same field drift.

**Fix:** [`join` on `meta.id` only, then `multiMap` back to the two channels BCFTOOLS_VIEW
expects](../../subworkflows/local/pharmcat_vcf_processing.nf). Standard nf-core pattern,
works regardless of other meta-field drift.

## Running the demo

```bash
bash tests/mwe/run_demo.sh
```

The script:

1. Runs both MWEs with the current working tree (the branch's fix) and
   asserts the fix works.
2. `git stash`-es the fix from the two patched files and re-runs to
   demonstrate the pre-fix bugs reproduce.
3. Restores the fix (also on Ctrl-C / errors via an `EXIT` trap).

Expected output ends with `All N/N assertions passed.`

## Requirements

- `git`, `bash`
- [`pixi`](https://pixi.sh) — everything else (`nextflow`, `samtools`,
  `htslib`, `bcftools`) is fetched on demand via `pixi exec`.

The harnesses run in Nextflow stub mode against the actual project
modules and subworkflow — no external test datasets, no real PharmCAT
container, ~30 seconds end-to-end.

## What each harness exercises

| MWE | What it forces |
|-----|----------------|
| `test_26_stage.nf` | Inspects the staged inputs in the `PHARMCAT_VCFPREPROCESSOR` workdir. Asserts 0 symlinks under the fix. |
| `test_27_drift.nf` | Two samples, VCF meta has extra `data_type` field, BED meta has only `id`, BED emission order *reversed*. Asserts each sample's BCFTOOLS_VIEW task receives its own BED under the fix. |

The dummy input files (`tests/mwe/inputs/`) are empty stubs created on
the fly — Nextflow's stub-mode tasks don't read their contents.
