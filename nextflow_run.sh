#!/bin/bash
#SBATCH --account=p212_tsd
#SBATCH --job-name=nextflow
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --time=2:00:00
#SBATCH --output=nextflow.log
#SBATCH --error=nextflow.err

module load Nextflow/24.04.2
nextflow run main.nf -profile p212 -resume
