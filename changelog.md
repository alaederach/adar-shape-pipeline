# Changelog

Records every change to `main.nf`, `step02.nf`, or a module, grouped by file, as a
`**Fix**`/`**Change**` pair explaining the failure mode and the resolution. Also
documents hard-won Nextflow gotchas — read this before touching the cache, output
globs, or join keys.

## Nextflow gotchas

- **Editing a `bin/` script does NOT invalidate the `-resume` cache.** Nextflow keys a
  task's cache on the process's inline `script:` text plus its declared inputs — *not*
  the contents of any external script the body calls (e.g.
  `python3 ${projectDir}/bin/plot_skyline.py`). So after editing `bin/plot_skyline.py`
  and re-running with `-resume`, `PLOT_SKYLINE` came back `Cached` and the published PNG
  was byte-identical — the edit silently had no effect. Symptoms: you change a `bin/`
  script, re-run with `-resume`, and the output doesn't change.
  **Workarounds:** (a) re-run without `-resume` (re-runs everything — slow); (b) delete
  just that task's work dir (`rm -rf work/<hash-prefix>*`) then `-resume` to force only
  that process to re-run; (c) for plot-only iteration, run the `bin/` script directly
  against the already-published profiles, bypassing Nextflow. A fresh run (no `-resume`)
  is never affected. To make this airtight you could add the script as a declared `path`
  input so its hash is tracked, but that complicates every call site — not currently done.

## 2026-06-23 — distribution + code hardening

### lib/Validate.groovy (new)
- **Change** — Shared fail-fast validation (auto-loaded from `lib/`). `Validate.params`
  checks refbase/readbase (valid, distinct nucleotides), `min_depth` (positive int),
  `min_mutation_separation` (≥0). `Validate.step01Row` / `step02Row` check required
  columns, `sample_type`/`library_type` against allowed sets, `min_rate` numeric in (0,1],
  candidates as int list, and that `targeted-amplicon` rows carry a primers file. Errors
  name the offending row/column with the expected values.

### main.nf, step02.nf
- **Fix** — Malformed samplesheets failed late and cryptically (a typo'd `sample_type`
  produced an id like `STK4_AluSq2_null_rep1` and ran ShapeMapper with the wrong args).
  Now both call `Validate.params` + the per-row validator before building `meta`, open
  the samplesheet with `checkIfExists`, resolve `fasta` with `checkIfExists`, and verify
  the `mod_fastq`/`unt_fastq` folders and primers file exist — all with clear messages.

### nextflow.config
- **Fix** — The `longleaf` profile hardcoded a personal cluster path (`/work/users/.../`),
  so it leaked a personal path and worked for no one else. Renamed to a generic `slurm`
  profile driven by `params.shapemapper_bin` and `params.cluster_before_script` (each with
  an env-var fallback: `ADAR_SHAPEMAPPER_BIN`, `ADAR_CLUSTER_BEFORE_SCRIPT`); the UNC
  Longleaf values are kept only as a commented example.
- **Change** — Added a `manifest{}` block; `nextflowVersion = '>=25.04.0, <25.10.0'` makes
  Nextflow fail fast with a clear message on an unsupported (e.g. NF 26) engine.

### modules/revert_snv.nf
- **Fix** — The whole-reference fallback (no `alu_start`/`alu_end`) took the first `@SQ`
  `LN:` from the BAM header, silently reverting only the first contig's range on a
  multi-contig reference. Now errors with a clear message if the reference has ≠1 contig
  and no explicit coordinates were given.

### bin/query_sorted_bam.py
- **Change** — `choose_best_threshold` now emits a stderr WARNING when even the loosest
  (≥1 edited position) threshold cannot reach `min_depth` reads in both BAMs — i.e. the
  cutoff is too stringent / data too shallow and the edited set (and its re-ShapeMapper
  profile) will be sparse. Previously it silently returned threshold 1.

## 2026-06-22 — step 02 deconvolution + skyline control plot

### step02.nf
- **Change** — Extended step 02 beyond read sorting to the full deconvolution: sort the
  edited (A→G) read set → revert the edits → re-run ShapeMapper on the reverted reads →
  normalize against the standard-SHAPE reference. New stages wire
  `REVERT_SNV` → `BAM_TO_FASTQ` → `SHAPEMAPPER` (alias `RESHAPE`, single-end) →
  `NORMALIZE_PROFILES` → `PLOT_SKYLINE`. Endpoint is the deconvolved + normalized
  reactivity profile of the edited molecules.
- **Change** — The re-ShapeMapper (`RESHAPE`) meta carries `min_rate`/`candidates` so the
  skyline plot can draw the cutoff line and sorted positions; `dp` (ADAR-DMSO) and `sp`
  (standard-SHAPE) profiles are joined back from `samples` by `ckey` for the plot.

### nextflow.config
- **Fix** — `step02.edited_min_depth` (read floor for threshold selection) was decoupled
  from `min_depth` (ShapeMapper's per-position HQ floor): 1000 vs 5000. The sorter would
  pick the highest "n-of-N edited" threshold meeting the *read* floor (e.g. ≥3-of-7,
  ~1600 reads), but the re-ShapeMapper then masked every position below the *depth* floor
  → the normalized profile came out all-`nan`. Resolution: coupled
  `edited_min_depth = unedited_min_depth = 5000` (= `min_depth`), forcing a threshold
  whose read count clears the HQ depth floor. Raise both together for deeper data.

### modules/normalize_profiles.nf
- **Change** — Also emit the normalized *reference* (`--normout`), not just the scaled
  edited profile, so the skyline plot compares two tracks produced on one normalization
  scale (matched `--normout`/`--scaleout` from a single normalize call).

### modules/plot_skyline.nf, bin/plot_skyline.py (new)
- **Change** — Two-panel skyline control figure (RNAvigate `plot_profile_skyline` style:
  `ax.plot(drawstyle="steps-mid")` + `fill_between(step="mid")`, matplotlib only).
  Top: bulk SHAPE (blue) vs edited SHAPE-ADAR (orange) reactivity. Bottom: ADAR-DMSO
  editing (red) vs non-ADAR DMSO control (grey, = standard-SHAPE `Untreated_rate`), with
  the sorting cutoff as a horizontal line and every sorted position as a vertical line
  spanning both panels. Sorted positions come from `--sort-positions` (candidates mode)
  or are derived as refbase positions ≥ `--cutoff` (min_rate mode).

### modules/shapemapper.nf
- **Change** — Added single-end (`--U`) input mode (`meta.fastq_mode == 'single_end'`)
  for the reverted-read re-ShapeMapper; folder mode (raw step-01 FASTQs) is unchanged.

### data/AluSq2 (bundled demo)
- **Fix** — Re-downsampled the demo data 100k → ~300k reads/file (~178 MB). At 100k the
  edited subset was too shallow to clear `min_depth=5000` even after the threshold fix, so
  the demo's normalized profile was empty. At ~300k the ≥2-of-7 edited subset reaches
  ~6.9k depth and produces a real profile at default settings (no per-run flags).
