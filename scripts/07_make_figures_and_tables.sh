#!/usr/bin/env bash
# Aggregate both cohorts' data/public/ results into the final manuscript figures/tables.
# Run locally (not on TSD) after downloading data/public/ for both p212 and p2521.
# Usage: scripts/07_make_figures_and_tables.sh [--force]
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"
cd "$REPO_ROOT"

[[ "${1:-}" == "--force" ]] && FORCE=1

outputs=(manuscript/figures_and_tables/km_curve.pdf manuscript/figures_and_tables/table_coxph.tsv)
should_run "${outputs[@]}" || exit 0

mkdir -p figures_and_tables
run_rscript src/make_figures_and_tables.R
require_files figures_and_tables/km_curve.pdf figures_and_tables/table_coxph.tsv

mkdir -p manuscript/figures_and_tables
mv figures_and_tables/* manuscript/figures_and_tables/
rmdir figures_and_tables

echo "Done. Outputs written to manuscript/figures_and_tables/"
echo "Render the manuscript manually with: cd manuscript && quarto render manuscript.qmd --to pdf"
