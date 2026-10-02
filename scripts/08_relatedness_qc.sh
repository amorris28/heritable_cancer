#!/usr/bin/env bash
#SBATCH --job-name=kingrelatedness
#SBATCH --account=p2521
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3G
#SBATCH --time=01:00:00
#SBATCH --output=kingrelatedness.log
set -euo pipefail
# sbatch copies this script into a spool dir, so BASH_SOURCE can't locate the repo; use SLURM_SUBMIT_DIR instead.
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  SCRIPT_DIR="${SLURM_SUBMIT_DIR}/scripts"
else
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
source "${SCRIPT_DIR}/lib/common.sh"

COHORT="${1:?Usage: $0 <p212|p2521> [--force]}"
[[ "${2:-}" == "--force" ]] && FORCE=1
load_cohort_config "$COHORT"
cd "$REPO_ROOT"

OUT="data/interim/king"
RELATED="$OUT/${CLUSTER}related_deg3"
RELATED_FILES=("$RELATED.kin" "$RELATED.kin0")
PLINK="../bin/plink"
KING="../bin/king"

if ! should_run "${RELATED_FILES[@]}"; then
  exit 0
fi

mkdir -p data/interim/king

STEP5="$OUT/$COHORT.pca.step5"
if [[ "$FORCE" != "1" && -f "$STEP5.bed" && -f "$STEP5.bim" && -f "$STEP5.fam" ]]; then
  echo "  [skip] PLINK preprocessing outputs already present: $STEP5.{bed,bim,fam}"
else
  [[ -f "$BED_FILE.bed" ]] || { echo "Missing $BED_FILE.bed" >&2; exit 1; }

  $PLINK --bfile "$BED_FILE" --geno 0.05 --make-bed --out "$OUT/$COHORT.pca.step1"
  $PLINK --bfile "$OUT/$COHORT.pca.step1" --chr 1-22 --indep-pairwise 1000 50 0.05 --make-bed --out "$OUT/$COHORT.pca.step2"
  $PLINK --bfile "$OUT/$COHORT.pca.step1" --extract "$OUT/$COHORT.pca.step2.prune.in" --make-bed --out "$OUT/$COHORT.pca.step3"
  $PLINK --bfile "$OUT/$COHORT.pca.step3" --make-set data/external/high-ld.txt --write-set --out "$OUT/$COHORT.pca.step4"
  $PLINK --bfile "$OUT/$COHORT.pca.step3" --exclude "$OUT/$COHORT.pca.step4.set" --make-bed --out "$STEP5"
fi
require_files "$STEP5.bed" "$STEP5.bim" "$STEP5.fam"

# KING defaults to half of the node's (logical) cores for --cpus, which can
# wildly oversubscribe an interactive session or a Slurm job that was only
# allocated 1-2 CPUs, causing massive thread contention instead of speedup.
# Pin it to what's actually available (falling back to 1 if undetectable).
KING_CPUS="${SLURM_CPUS_PER_TASK:-$(nproc 2>/dev/null || echo 1)}"

$KING -b "$OUT/$COHORT.pca.step5.bed" \
  --related --degree 3 --cpus "$KING_CPUS" --prefix "$RELATED"

require_files "${RELATED_FILES[@]}"

echo "Done. Output: $RELATED"
