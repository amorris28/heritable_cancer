# Analysis pipeline

Each script is safe to (re)run on its own: it checks whether its output already
exists in `data/interim/`, `data/processed/`, or `data/public/` and skips work
unless `FORCE=1` is set (or `--force` passed). Nothing reruns automatically —
you decide what to rerun, so editing an early step never silently retriggers
an expensive downstream step (e.g. the bootstrap in `05_fit_models.sh`).

## Order of operations, per cohort (`p212` = NFCS/MoBA, `p2521` = HUSK)

1. `sbatch scripts/01_prepare_genotypes_<cohort>.sh` — one-time, long-running
   (~4h) genotype prep on TSD: VCF/BIM harmonization against the PHS601 model,
   SNP filtering, merging, and PHS scoring. Submitted via `sbatch` because it's
   slow and only needs to run once (or after a genotype/model update).
2. `scripts/02_preprocess_metadata.sh <cohort>` — builds `data/processed/phenotypes.tsv`
   from cohort-specific phenotype/registry sources + the PHS score. Fast, safe
   to rerun whenever phenotype logic changes.
3. `scripts/03_p212_relatedness.sh` (p212 only) — build KING retained-ID lists
   from the external kinship matrix. For p2521 the lists come from the external
   run; confirm they are available and built on the same genotype release.
4. `sbatch --account=<cohort> scripts/04_pca.sh <cohort>` — LD-pruned FlashPCA
   genetic PCs. Steps 3 and 4 are independent of each other and of step 2, but
   all must finish before model fitting.
5. `scripts/05_fit_models.sh <cohort>` — primary Cox models (excluding
   duplicates/MZ twins) + bootstrap performance metrics (999 iterations by
   default; override with `BOOTSTRAP_ITERS=99` for quick dev iteration), age-band
   HRs, and the genetic-QC sensitivity analyses.
6. `scripts/06_plot_figures.sh <cohort>` — KM/absolute-risk curve objects.
7. `scripts/07_incidence_analysis.sh <cohort>` — population incidence-rate comparison.
8. `scripts/08_wgs_analysis.sh` — rare-variant (WGS) analysis, p212 only.
9. `scripts/09_phs_portability.sh <cohort>` — PHS portability in the original
   score units (not used by the manuscript tables/figures).
10. `scripts/10_calibration_frozen.sh <cohort>` — calibration with the published
    PHS601 weights (not used by the manuscript tables/figures).
11. `scripts/11_calibration.sh <cohort>` — internal absolute-risk calibration
    (supplementary figure).
12. `scripts/12_make_figures_and_tables.sh` — run **locally**, after both cohorts'
    `data/public/` results are available, to build `manuscript/figures_and_tables/`.

After steps 2 and 3, run `bash scripts/13_sample_counts.sh <p212|p2521>` to add
the age-40+ before/after KING counts to the report created during preprocessing.
This is independent of steps 4-11 and can be run before model fitting. It only
reads the existing phenotype file and KING retained-ID lists; it does not rerun
preprocessing, fit models, or run bootstraps. If preprocessing logic has changed
or the report is missing, run `scripts/02_preprocess_metadata.sh <cohort>
--force` first. The output is `data/public/<cohort>_sample_counts.tsv`.

Manuscript rendering itself stays manual and separate:
```
quarto render manuscript/manuscript.qmd --to pdf \
  --metadata csl=journal-of-clinical-oncology.csl --execute-dir "$(pwd)"
```
`--execute-dir` must point at the repo root (run this from there) — the R
chunks in `manuscript.qmd` read `src/` and `data/public/` relative to it, but
Quarto otherwise executes code relative to the `.qmd`'s own directory.

## Relatedness and PC sensitivity analyses

The primary Cox models in `scripts/05_fit_models.sh <cohort>` use all eligible
subjects except one member of each KING 0-degree (duplicate/MZ) pair, with no
genetic PCs. Sensitivity
analyses estimate the PHS coefficient for each of three endpoints (any PCa,
metastatic PCa and PCa-specific death) across a 3-by-3 grid:

- retained set: KING 0-degree (primary), 2-degree and 3-degree
  retained sets;
- PC adjustment: none, first 5 PCs, or first 10 PCs.

Results, including kinship cutoff, PC count, sample size, event count, hazard
ratio and confidence interval, are written to
`data/public/<cohort>_genetic_qc_sensitivity.tsv`. The 0-degree list uses a
kinship cutoff of 0.35333 (duplicates/monozygotic-level); 2-degree and
3-degree lists use cutoffs 0.0883883 and 0.0441942, respectively.

- p2521: genotypes are the v260302 imputation (`config/p2521.env`); KING
  retained-ID lists must come from the same release. They come from the external run,
  `/ess/p2521/cluster/users/espehage/BiobankNorge_AD/results/husk/all_chr.cutoff_{0deg,2deg,3deg}.king.cutoff.in.id`.
- p212: run `scripts/03_p212_relatedness.sh` first. It converts the
  colleague's KING text matrix to `.king.bin` (`scripts/lib/king_text_to_bin.py`)
  and runs `plink2 --king-cutoff` at each threshold, writing
  `data/interim/relatedness/p212_{0deg,2deg,3deg}.king.cutoff.in.id`. Override
  inputs with `KING_TEXT`, `KING_IDS` or `KING_BFILE`.

The retained-ID lists provide progressively less-related subsets, not
family/cluster identifiers, so they do not support a full-cohort
relatedness-adjusted Cox model. PC scores must cover every analysis sample;
the model script checks the IDs and requires at least 10 PC columns.

## Participant flow counts

`data/public/<cohort>_sample_counts.tsv` is created by step 2 and records the
source cohort, male participants, overlap with the post-QC genotype and
PHS-scoring files, complete phenotype/PHS sets, and the age-40+ primary-model
input. Step 13 adds the exact before/after KING counts, including the number
removed by the 0-degree duplicate/MZ filter. For NFCS the preprocessing report
also includes WGS-flagged
participants in the source cohort, male WGS participants, and the WGS analysis
set. These are distinct denominators; use the stage definitions in the file
when reporting them. The previously cited 2,576 sequenced participants is an
overall sequencing-cohort count and should not be presented as the analytic
male WGS sample size.

## PCA

Submit `scripts/04_pca.sh <p212|p2521>` with the matching Slurm account, for
example:

```sh
sbatch --account=p2521 scripts/04_pca.sh p2521
sbatch --account=p212 scripts/04_pca.sh p212
```

The script uses the configured cohort genotype BED set, excludes build-matched
high-LD intervals during `--indep-pairwise 1000 50 0.05` pruning, and runs
FlashPCA using the local `../bin/plink` and `../bin/flashpca_x86-64` binaries.
It passes the build-matched PLINK range file directly to
`--exclude range`: `scripts/resources/high-LD-regions-hg38-GRCh38.txt` for
GRCh38 (from
[plinkQC](https://github.com/cran/plinkQC/blob/master/inst/extdata/high-LD-regions-hg38-GRCh38.txt))
and `scripts/resources/high-LD-regions-hg19-GRCh37.txt` for GRCh37. The files
are used as-is without coordinate or chromosome-name conversion.

PLINK intermediates and FlashPCA outputs are written under
`data/interim/pca/`. The `pcs_<cohort>.txt` file contains the sample PC scores;
loadings and means/standard deviations are also retained for later projection.

## Config

Per-cohort paths/params (TSD genotype file locations, container images, SLURM
account) live in `config/p212.env` and `config/p2521.env` and are sourced
automatically by the scripts — nothing to edit unless those TSD paths change.

## Forcing a rerun

Every script accepts `--force` (or `FORCE=1` env var) to ignore existing
outputs and regenerate them. Otherwise it prints `[skip]` and exits.
