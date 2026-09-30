#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

COHORT="${1:?Usage: $0 <p212|p2521> [--force]}"
[[ "${2:-}" == "--force" ]] && FORCE=1
load_cohort_config "$COHORT"
cd "$REPO_ROOT"

OUT="data/interim/king"
RELATED="$OUT/${CLUSTER}related_deg3"
RELATED_FILES=("$RELATED.kin" "$RELATED.kin0")
PLINK="../bin/plink"
PLINK="../bin/king"

if ! should_run "${RELATED_FILES[@]}"; then
  exit 0
fi

[[ -f "$BED_FILE" ]] || { echo "Missing $BED_FILE" >&2; exit 1; }

mkdir -p data/interim/king

$PLINK --bfile $BED_FILE --geno 0.05 --make-bed --out $OUT/$COHORT.pca.step1
$PLINK --bfile $OUT/$COHORT.pca.step1 --chr 1-22 --indep-pairwise 1000 50 0.05 --make-bed --out $OUT/$COHORT.pca.step2
$PLINK --bfile $OUT/$COHORT.pca.step1 --extract $OUT/$COHORT.pca.step2.prune.in --make-bed --out $OUT/$COHORT.pca.step3
$PLINK --bfile $OUT/$COHORT.pca.step3 --make-set $OUT/high-ld.txt --write-set --out $OUT/$COHORT.pca.step4
$PLINK --bfile $OUT/$COHORT.pca.step3 --exclude $OUT/$COHORT.pca.step4.set --make-bed --out $OUT/$COHORT.pca.step5

$KING -b "$$OUT/$COHORT.pca.step5.bed" \
  --related --degree 3 --prefix "$RELATED"

require_files "${RELATED_FILES[@]}"

echo "Done. Output: $RELATED"
