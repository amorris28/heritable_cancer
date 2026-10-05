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

## Relatedness QC

Run `sbatch scripts/08_relatedness_qc.sh <p212|p2521>` to identify close
relatives with KING. Its KING input is restricted to autosomes and variants
with at most 5% missingness, but is not LD-pruned and does not exclude
high-LD regions. This keeps the relatedness analysis separate from the
LD-pruned genotype set that will be prepared for PCA. Results from the earlier
pruned workflow use a different output prefix and are not reused.

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
