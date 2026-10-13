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

## p212 KING relatedness

Submit the KING jobs from the repository root, starting with
`king/make_king.job`. After it completes, submit
`king/king_cutoff_{0,1,2,3}.job` to generate the unrelated-sample lists at
kinship cutoffs 0.353553, 0.176777, 0.0883883, and 0.0441942, respectively.
The default full-genome PLINK binary prefix is
`/ess/p212/cluster/moba_qc_imputation/BA/release6/moba_imputed_all`; override
it with `KING_BFILE` when submitting if needed. Results are written to
`data/interim/relatedness/p212/`.

## p2521 Relatedness sensitivity analysis

The p2521 KING cutoff lists are supplied by the external relatedness run.
`scripts/03_fit_models.sh p2521` uses the 0-degree retained list (kinship
cutoff 0.35333) as its primary analysis set, removing potential
duplicates/monozygotic-level pairs while retaining one sample from each pair.
It also estimates the PHS Cox coefficients in the 0-degree, 2-degree and
3-degree retained sets (kinship cutoffs 0.35333, 0.0883883 and 0.0441942). The
estimates are written to `data/public/p2521_relatedness_sensitivity.tsv`.

These ID lists provide unrelated subsets, not family/cluster identifiers, so
they do not support a full-cohort relatedness-adjusted Cox model. The
sensitivity models should be interpreted as robustness checks against
excluding increasingly distant relatives. The required input lists are
`all_chr.cutoff_{0deg,2deg,3deg}.king.cutoff.in.id` under
`/ess/p2521/cluster/users/espenhage/BioBankNorge_AD/results/husk`.

## PCA

Submit `scripts/09_pca.sh <p212|p2521>` with the matching Slurm account, for
example:

```sh
sbatch --account=p2521 scripts/09_pca.sh p2521
sbatch --account=p212 scripts/09_pca.sh p212
```

The script uses the autosomal KING input when available (otherwise it creates
one with `--geno 0.05`), excludes build-matched high-LD intervals during
`--indep-pairwise 1000 50 0.05` pruning, and runs FlashPCA using the local
`../bin/plink` and `../bin/flashpca_x86-64` binaries. It passes the build-matched PLINK range file directly to
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
