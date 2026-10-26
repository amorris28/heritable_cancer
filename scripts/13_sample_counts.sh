#!/usr/bin/env bash
# Add KING retained/excluded counts to the per-cohort preprocessing report.
# Usage: bash scripts/13_sample_counts.sh <p212|p2521>
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"
COHORT="${1:?Usage: $0 <p212|p2521>}"
load_cohort_config "$COHORT"
cd "$REPO_ROOT"

if [[ "$COHORT" == "p2521" ]]; then
  cohort_label=p2521
  relatedness_dir=/ess/p2521/cluster/users/espehage/BiobankNorge_AD/results/husk
  relatedness_files=(
    "$relatedness_dir/all_chr.cutoff_0deg.king.cutoff.in.id"
    "$relatedness_dir/all_chr.cutoff_2deg.king.cutoff.in.id"
    "$relatedness_dir/all_chr.cutoff_3deg.king.cutoff.in.id"
  )
else
  cohort_label=p212
  relatedness_files=(
    "data/interim/relatedness/p212_0deg.king.cutoff.in.id"
    "data/interim/relatedness/p212_2deg.king.cutoff.in.id"
    "data/interim/relatedness/p212_3deg.king.cutoff.in.id"
  )
fi

[[ -f data/processed/phenotypes.tsv ]] || {
  echo "Missing data/processed/phenotypes.tsv — run scripts/02_preprocess_metadata.sh $COHORT first" >&2
  exit 1
}
[[ -f "data/public/${CLUSTER}_sample_counts.tsv" ]] || {
  echo "Missing data/public/${CLUSTER}_sample_counts.tsv — rerun scripts/02_preprocess_metadata.sh $COHORT --force first" >&2
  exit 1
}
for file in "${relatedness_files[@]}"; do
  [[ -f "$file" ]] || { echo "Missing relatedness keep list: $file" >&2; exit 1; }
done

run_rscript src/sample-counts.R \
  data/processed/phenotypes.tsv \
  "$cohort_label" \
  "${relatedness_files[@]}" \
  "data/public/${CLUSTER}_sample_counts.tsv"

echo "Done. Added KING counts to data/public/${CLUSTER}_sample_counts.tsv"
