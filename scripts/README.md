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

## Relatedness sensitivity analysis

`scripts/03_fit_models.sh <cohort>` uses KING cutoff retained-ID lists
(`plink2 --king-cutoff` `.in.id` files). The 0-degree list (kinship cutoff
0.35333) defines the primary analysis set, removing duplicate/monozygotic-level
pairs while retaining one sample from each pair. The PHS Cox coefficients are
also estimated in the 0-degree, 2-degree and 3-degree retained sets (cutoffs
0.35333, 0.0883883 and 0.0441942) and written to
`data/public/<cohort>_relatedness_sensitivity.tsv`.

- p2521: lists come from the external run,
  `/ess/p2521/cluster/users/espehage/BiobankNorge_AD/results/husk/all_chr.cutoff_{0deg,2deg,3deg}.king.cutoff.in.id`.
- p212: run `scripts/08_p212_relatedness.sh` first. It converts the
  colleague's KING text matrix to `.king.bin` (`scripts/lib/king_text_to_bin.py`)
  and runs `plink2 --king-cutoff` at each threshold, writing
  `data/interim/relatedness/p212_{0deg,2deg,3deg}.king.cutoff.in.id`. Override
  inputs with `KING_TEXT`, `KING_IDS` or `KING_BFILE`.

These ID lists provide unrelated subsets, not family/cluster identifiers, so
they do not support a full-cohort relatedness-adjusted Cox model. The
sensitivity models are robustness checks against excluding increasingly
distant relatives.

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
