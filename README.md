# Prostate Cancer Polygenic Hazard Score in Norwegian Cohorts

## Data Availability

Polygenic hazard score (PHS) model parameters for PHS601 were acquired
through personal communication with the authors.

Genotype and phenotype data were accessed through the University
of Oslo's Tjenester for Sensitiv Data (TSD).

## Required software

Analyses are run as a series of plain bash/R/python scripts. 
See `scripts/README.md` for the full list and order of operations.
`bash`, `R` (>= 4.4), `python3`, `plink`/`plink2`, and `quarto` (for the
manuscript) are required, along with `singularity` on TSD.

All R, python3, and plink calls were performed using
[the COMORMENT containers](https://github.com/comorment/containers) with either
singularity or docker. Analyses were developed using the COMORMENT containers
[version 1.9.0](https://github.com/comorment/containers/releases/tag/v1.9.0). Compiling the manuscript requires additional packages included in the provided Dockerfile.

## Commands to reproduce analysis

Run the numbered scripts in `scripts/` in order, once per TSD cohort project
(`p212` and `p2521`). Each script skips its work if its
output already exists, so reruns only recompute what you ask for:

```
sbatch scripts/01_prepare_genotypes_p212.sh   # one-time, long-running genotype prep (per cohort)
scripts/02_preprocess_metadata.sh p212        # phenotype preprocessing
scripts/03_fit_models.sh p212                 # Cox models + bootstrap performance
scripts/04_plot_figures.sh p212               # KM/absolute-risk curves
scripts/05_incidence_analysis.sh p212         # population incidence comparison
scripts/06_wgs_analysis.sh                    # rare-variant (WGS) analysis, p212 only
```

Repeat steps 2-5 for `p2521`. Once `data/public/` has been exported from both
TSD projects to a local machine, build the manuscript figures/tables and
render the manuscript:

```
scripts/07_make_figures_and_tables.sh
cd manuscript && quarto render manuscript.qmd --to pdf
```

See `scripts/README.md` for details, cohort-specific config (`config/p212.env`,
`config/p2521.env`), and how to force a step to rerun.

## Project Organization

```
Project/
│                         
├── scripts/           <- Numbered, idempotent scripts that run the analysis
│                         (see scripts/README.md for order of operations)
│                         
├── config/            <- Per-cohort (p212, p2521) paths/params, sourced by scripts/
│                         
├── Makefile           <- Makefile with convenience commands 
│
├── README.md          <- The top-level README
│                         
├── data
│   ├── external       <- Data from third-party sources
│   ├── raw            <- Raw data
│   ├── interim        <- Intermediate data that has been transformed
│   ├── processed      <- The final, canonical data sets for modeling
│   └── public         <- Summary data that can be download from TSD
│
├── manuscript         <- Manuscript documents
│   ├── manuscript.qmd <- Quarto markdown document containing the main text
│   ├── reference.docx <- Reference file for formatting Quarto docx output
│   └── references.bib <- Bibtex references for citing in the text
│
├── src                <- Source code for R, Python, and bash scripts used 
│                         in this project (wrapped by scripts/)
│                         
└── archive            <- Retired Nextflow pipeline and superseded analyses,
                          kept for reference only
```
