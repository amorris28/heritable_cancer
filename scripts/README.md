# Analysis pipeline

Each script is safe to (re)run on its own: it checks whether its output already
exists in `data/interim/`, `data/processed/`, or `data/public/` and skips work
unless `FORCE=1` is set (or `--force` passed). Nothing reruns automatically —
you decide what to rerun, so editing an early step never silently retriggers
an expensive downstream step (e.g. the bootstrap in `03_fit_models.sh`).

## Order of operations, per cohort (`p212` = NFCS/MoBA, `p2521` = HUSK)

1. `sbatch scripts/01_prepare_genotypes_<cohort>.sh` — one-time, long-running
   (~4h) genotype prep on TSD: VCF/BIM harmonization against the PHS601 model,
   SNP filtering, merging, and PHS scoring. Submitted via `sbatch` because it's
   slow and only needs to run once (or after a genotype/model update).
2. `scripts/02_preprocess_metadata.sh <cohort>` — builds `data/processed/phenotypes.tsv`
   from cohort-specific phenotype/registry sources + the PHS score. Fast, safe
   to rerun whenever phenotype logic changes.
3. `scripts/03_fit_models.sh <cohort>` — Cox models + bootstrap performance
   metrics (999 iterations by default; override with `BOOTSTRAP_ITERS=99` for
   quick dev iteration).
4. `scripts/04_plot_figures.sh <cohort>` — KM/absolute-risk curve objects.
5. `scripts/05_incidence_analysis.sh <cohort>` — population incidence-rate comparison.
6. `scripts/06_wgs_analysis.sh` — rare-variant (WGS) analysis, p212 only.
7. `scripts/07_make_figures_and_tables.sh` — run **locally**, after both cohorts'
   `data/public/` results are available, to build `manuscript/figures_and_tables/`.

Manuscript rendering itself stays manual and separate:
```
quarto render manuscript/manuscript.qmd --to pdf \
  --metadata csl=journal-of-clinical-oncology.csl --execute-dir "$(pwd)"
```
`--execute-dir` must point at the repo root (run this from there) — the R
chunks in `manuscript.qmd` read `src/` and `data/public/` relative to it, but
Quarto otherwise executes code relative to the `.qmd`'s own directory.

## Config

Per-cohort paths/params (TSD genotype file locations, container images, SLURM
account) live in `config/p212.env` and `config/p2521.env` and are sourced
automatically by the scripts — nothing to edit unless those TSD paths change.

## Forcing a rerun

Every script accepts `--force` (or `FORCE=1` env var) to ignore existing
outputs and regenerate them. Otherwise it prints `[skip]` and exits.
