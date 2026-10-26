#!/usr/bin/env bash
# Phenotype preprocessing — fast, safe to rerun often when you tweak phenotype logic.
# Usage: scripts/02_preprocess_metadata.sh <p212|p2521> [--force]
# On TSD, run `module load R/4.4.2-gfbf-2024a` first (or add it to your shell rc).
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

COHORT="${1:?Usage: $0 <p212|p2521> [--force]}"
[[ "${2:-}" == "--force" ]] && FORCE=1
load_cohort_config "$COHORT"
cd "$REPO_ROOT"

PHENOTYPES=data/processed/phenotypes.tsv
SCORE_FILE=data/processed/${CLUSTER}_phs601.profile
PHS_BFILE=data/interim/merge_bfiles/${CLUSTER}_phs601
SAMPLE_COUNTS="data/public/${CLUSTER}_sample_counts.tsv"
AGE_SUMMARY="data/public/${CLUSTER}_age_summary.tsv"
KING_0DEG="$(king_keep_file "$CLUSTER" 0deg)"

if ! should_run "$PHENOTYPES" "$SAMPLE_COUNTS" "$AGE_SUMMARY"; then
  exit 0
fi

[[ -f "$KING_0DEG" ]] || { echo "Missing $KING_0DEG (needed for Table 1 summaries) — run scripts/03_p212_relatedness.sh first (p212) or check the HUSK KING path" >&2; exit 1; }
[[ -f "$SCORE_FILE" ]] || { echo "Missing $SCORE_FILE — run scripts/01_prepare_genotypes_${CLUSTER}.sh first" >&2; exit 1; }
if [[ "$HAS_HOXB13" == "1" ]]; then
  [[ -f data/processed/tmp_rs138213197.raw ]] || { echo "Missing data/processed/tmp_rs138213197.raw — run scripts/01_prepare_genotypes_${CLUSTER}.sh first" >&2; exit 1; }
  cp data/processed/tmp_rs138213197.raw tmp_rs138213197.raw
fi

run_rscript "$PREPROCESS_SCRIPT" "$SCORE_FILE" FALSE "$BED_FILE" "$PHS_BFILE" "$SAMPLE_COUNTS" "$KING_0DEG"
require_files phenotypes.tsv age_summary.tsv "$SAMPLE_COUNTS"
if [[ "$HAS_HOXB13" == "1" ]]; then
  require_files all_variants.tsv age_summary_wgs.tsv
fi

mv phenotypes.tsv "$PHENOTYPES"
mv age_summary.tsv "data/public/${CLUSTER}_age_summary.tsv"
if [[ "$HAS_HOXB13" == "1" ]]; then
  mv all_variants.tsv data/public/all_variants.tsv
  mv age_summary_wgs.tsv "data/public/${CLUSTER}_age_summary_wgs.tsv"
fi
rm -f tmp_rs138213197.raw

echo "Done. Output: $PHENOTYPES"
