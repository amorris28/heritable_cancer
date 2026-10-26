#!/usr/bin/env bash
# Add KING retained/excluded counts to the per-cohort preprocessing report.
# Usage: bash scripts/13_sample_counts.sh <p212|p2521>
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"
COHORT="${1:?Usage: $0 <p212|p2521>}"
load_cohort_config "$COHORT"
cd "$REPO_ROOT"

cohort_label="$COHORT"
relatedness_files=(
  "$(king_keep_file "$COHORT" 0deg)"
  "$(king_keep_file "$COHORT" 2deg)"
  "$(king_keep_file "$COHORT" 3deg)"
)

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
