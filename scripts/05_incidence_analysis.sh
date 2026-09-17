#!/usr/bin/env bash
# Population incidence-rate comparison (Norwegian Cancer Registry ASIRs vs cohort).
# Usage: scripts/05_incidence_analysis.sh <p212|p2521> [--force]
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
  "data/public/${CLUSTER}_cumulative_comparison.rds"
  "data/public/${CLUSTER}_stratified_incidence.rds"
  "data/public/${CLUSTER}_cumulative_comparison_metastatic.rds"
  "data/public/${CLUSTER}_stratified_incidence_metastatic.rds"
)
should_run "${outputs[@]}" || exit 0

Rscript src/fit-models-plot-figures-incidence-rate.R "$PHENOTYPES" "$CLUSTER"

mv "${CLUSTER}_cumulative_comparison.rds" "${CLUSTER}_stratified_incidence.rds" \
   "${CLUSTER}_cumulative_comparison_metastatic.rds" "${CLUSTER}_stratified_incidence_metastatic.rds" \
   data/public/

echo "Done. Outputs: ${outputs[*]}"
