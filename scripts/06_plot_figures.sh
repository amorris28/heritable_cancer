#!/usr/bin/env bash
# KM/absolute-risk curve objects for the manuscript figures.
# Usage: scripts/04_plot_figures.sh <p212|p2521> [--force]
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

COHORT="${1:?Usage: $0 <p212|p2521> [--force]}"
[[ "${2:-}" == "--force" ]] && FORCE=1
load_cohort_config "$COHORT"
cd "$REPO_ROOT"

PHENOTYPES=data/processed/phenotypes.tsv
[[ -f "$PHENOTYPES" ]] || { echo "Missing $PHENOTYPES — run scripts/02_preprocess_metadata.sh $CLUSTER first" >&2; exit 1; }

outputs=(
  "data/public/${CLUSTER}_any_pca_abs_risk.rds"
  "data/public/${CLUSTER}_any_pca_risk_table.rds"
  "data/public/${CLUSTER}_metastatic_abs_risk.rds"
  "data/public/${CLUSTER}_metastatic_risk_table.rds"
  "data/public/${CLUSTER}_zph_tt.rds"
  "data/public/${CLUSTER}_zph_tt_aggr.rds"
  "data/public/${CLUSTER}_cuminc_at_age.tsv"
)
should_run "${outputs[@]}" || exit 0

run_rscript src/plot-figures.R "$PHENOTYPES" "$CLUSTER"
require_files "${CLUSTER}_any_pca_abs_risk.rds" "${CLUSTER}_any_pca_risk_table.rds" \
  "${CLUSTER}_metastatic_abs_risk.rds" "${CLUSTER}_metastatic_risk_table.rds" \
  "${CLUSTER}_zph_tt.rds" "${CLUSTER}_zph_tt_aggr.rds" "${CLUSTER}_cuminc_at_age.tsv"

mv "${CLUSTER}_any_pca_abs_risk.rds" "${CLUSTER}_any_pca_risk_table.rds" \
   "${CLUSTER}_metastatic_abs_risk.rds" "${CLUSTER}_metastatic_risk_table.rds" \
   "${CLUSTER}_zph_tt.rds" "${CLUSTER}_zph_tt_aggr.rds" "${CLUSTER}_cuminc_at_age.tsv" data/public/

echo "Done. Outputs: ${outputs[*]}"
