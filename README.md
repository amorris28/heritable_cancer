# Prostate Cancer Polygenic Hazard Score in Norwegian Cohorts

## Data Availability

Polygenic hazard score (PHS) model parameters for PHS601 were acquired
through personal communication with the authors.

Genotype and phenotype data were accessed through the University
of Oslo's Tjenester for Sensitiv Data (TSD).

## Required software

Recreating the analysis requires a recent version of [Nextflow](https://github.com/nextflow-io/nextflow).
Analyses were developed using Nextflow [version 24.04.2](https://github.com/nextflow-io/nextflow/releases/tag/v24.04.2).

All R, python3, and plink calls were performed using
[the COMORMENT containers](https://github.com/comorment/containers) with either 
singularity or docker. Analyses were developed using the COMORMENT containers 
[version 1.9.0](https://github.com/comorment/containers/releases/tag/v1.9.0).

## Commands to reproduce analysis

To run the analysis, run `main.nf` with nextflow. Passing a `-profile` is
_required_. Profiles are available for TSD projects `p2521` and `p212` as 
well as `local` for compiling of the manuscript. It is recommended to always
run with `-resume` so that cached files are re-used. Here is an example command:

```
nextflow run main.nf -profile p2521 -resume
```

The final output of `-profile local` (once all figures and tables have been
exported from TSD to local) is a copy of the compiled manuscript in `.docx`
and `.pdf` format.

## Project Organization

```
Project/
│                         
├── main.nf            <- Nextflow workflow file 
│                         
├── nextflow.config    <- Nextflow config with project-specific params
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
└── src                <- Source code for R, Python, and bash scripts used 
                          in this project
```

