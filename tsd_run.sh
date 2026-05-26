#!/bin/bash
# tsd_run.sh
# Generic job script executed on TSD via the HPC API.
# NOTE: SBATCH directives here are ignored — resource allocation is controlled
# by the submit_job() call in run_pipeline.py.
#
# Required variables (prepended by run_pipeline.py when generating per-project scripts):
#   PROJECT       e.g. p212
#   PROFILE       e.g. p212 or p2521
#   PIPELINE_DIR  e.g. /tsd/p212/home/p212-andrewmo/heritable_cancer
#   IMPORT_DIR    e.g. /tsd/p212/data/durable/file-import/p212-member-group
#   EXPORT_DIR    e.g. /tsd/p212/data/durable/file-export
set -euo pipefail

# Remove this uploaded job script from the import directory now that it has
# started. Safe on Linux — bash has already loaded it; the inode persists.
rm -f "$0"

mkdir -p "$PIPELINE_DIR"
cd "$PIPELINE_DIR"

echo "[$PROJECT] Unzipping pipeline..."
unzip -o "$IMPORT_DIR/pipeline.zip"

echo "[$PROJECT] Installing hazrd..."
module load R/4.4.2-gfbf-2024a
Rscript -e 'install.packages("hazrd_latest.tar.gz", repos=NULL, type="source", dependencies=FALSE)'

echo "[$PROJECT] Running Nextflow with profile: $PROFILE"
module load Nextflow/25.04.6
export NXF_HOME="$PIPELINE_DIR/.nextflow"
export NXF_WORK="$PIPELINE_DIR/work"
nextflow run main.nf -profile "$PROFILE" -resume

echo "[$PROJECT] Packaging results..."
zip -r public.zip data/public/*.tsv data/public/*.rds data/public/*.pdf data/public/*.png

echo "[$PROJECT] Exporting results..."
cp public.zip "$EXPORT_DIR/public_${PROJECT}.zip"

echo "[$PROJECT] Done."
