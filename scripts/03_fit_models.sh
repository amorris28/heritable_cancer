#!/usr/bin/env bash
# Cox models + bootstrap performance metrics — the slow interactive step (bootstrap_iters, default 999).
# Usage: scripts/03_fit_models.sh <p212|p2521> [--force]
# BOOTSTRAP_ITERS=99 scripts/03_fit_models.sh p212   # quick dev iteration with fewer bootstraps
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
  "data/public/${CLUSTER}_phs_zph_tidy.tsv"
  "data/public/${CLUSTER}_any_pca_performance.tsv"
  "data/public/${CLUSTER}_metastatic_performance.tsv"
  "data/public/${CLUSTER}_pca_death_performance.tsv"
  "data/public/${CLUSTER}_sample_size.tsv"
  "data/public/${CLUSTER}_coxph.tsv"
)

relatedness_files=()
if [[ "$COHORT" == "p2521" ]]; then
  relatedness_dir=/ess/p2521/cluster/users/espehage/BiobankNorge_AD/results/husk
  relatedness_files=(
    "$relatedness_dir/all_chr.cutoff_0deg.king.cutoff.in.id"
    "$relatedness_dir/all_chr.cutoff_2deg.king.cutoff.in.id"
    "$relatedness_dir/all_chr.cutoff_3deg.king.cutoff.in.id"
  )
  outputs+=("data/public/${CLUSTER}_relatedness_sensitivity.tsv")
  for f in "${relatedness_files[@]}"; do
    [[ -f "$f" ]] || { echo "Missing relatedness keep list: $f" >&2; exit 1; }
  done
fi

should_run "${outputs[@]}" || exit 0

r_args=("$PHENOTYPES" TRUE "$CLUSTER" "$BOOTSTRAP_ITERS")
if [[ "$COHORT" == "p2521" ]]; then
  r_args+=("${relatedness_files[@]}")
fi
run_rscript src/fit-models.R "${r_args[@]}"
require_files "${CLUSTER}_phs_zph_tidy.tsv" "${CLUSTER}_any_pca_performance.tsv" \
  "${CLUSTER}_metastatic_performance.tsv" "${CLUSTER}_pca_death_performance.tsv" \
  "${CLUSTER}_sample_size.tsv" "${CLUSTER}_coxph.tsv"

if [[ "$COHORT" == "p2521" ]]; then
  require_files "${CLUSTER}_relatedness_sensitivity.tsv"
  relatedness_output="${CLUSTER}_relatedness_sensitivity.tsv"
else
  relatedness_output=
fi

move_files=(
  "${CLUSTER}_phs_zph_tidy.tsv"
  "${CLUSTER}_any_pca_performance.tsv"
  "${CLUSTER}_metastatic_performance.tsv"
  "${CLUSTER}_pca_death_performance.tsv"
  "${CLUSTER}_sample_size.tsv"
  "${CLUSTER}_coxph.tsv"
)
if [[ -n "$relatedness_output" ]]; then
  move_files+=("$relatedness_output")
fi
mv "${move_files[@]}" data/public/

echo "Done. Outputs: ${outputs[*]}"
