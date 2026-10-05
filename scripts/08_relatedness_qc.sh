#!/usr/bin/env bash
#SBATCH --job-name=kingrelatedness
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --time=24:00:00
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
RELATED="$OUT/${CLUSTER}_related_deg3_unpruned"
RELATED_FILES=("$RELATED.kin" "$RELATED.kin0")
PLINK="../bin/plink"
KING="../bin/king"

if ! should_run "${RELATED_FILES[@]}"; then
  exit 0
fi

mkdir -p data/interim/king

KING_BFILE="$OUT/$COHORT.king_thin"
if [[ "$FORCE" != "1" && -f "$KING_BFILE.bed" && -f "$KING_BFILE.bim" && -f "$KING_BFILE.fam" ]]; then
  echo "  [skip] KING input files already present: $KING_BFILE.{bed,bim,fam}"
else
  [[ -f "$BED_FILE.bed" ]] || { echo "Missing $BED_FILE.bed" >&2; exit 1; }

  # KING cost scales with samples^2 x SNPs; ~6M imputed SNPs is far more than
  # needed (a few hundred thousand common SNPs resolve 3rd-degree kinship).
  # Random thinning (no LD pruning) keeps KING's assumptions intact.
  $PLINK --bfile "$BED_FILE" --chr 1-22 --geno 0.05 --maf 0.05 \
    --thin-count 400000 --seed 1 --make-bed --out "$KING_BFILE"
fi
require_files "$KING_BFILE.bed" "$KING_BFILE.bim" "$KING_BFILE.fam"

# KING defaults to half of the node's (logical) cores for --cpus, which can
# wildly oversubscribe an interactive session or a Slurm job that was only
# allocated 1-2 CPUs, causing massive thread contention instead of speedup.
# Pin it to what's actually available (falling back to 1 if undetectable).
KING_CPUS="${SLURM_CPUS_PER_TASK:-$(nproc 2>/dev/null || echo 1)}"

$KING -b "$KING_BFILE.bed" \
  --related --degree 3 --cpus "$KING_CPUS" --prefix "$RELATED"

require_files "${RELATED_FILES[@]}"

echo "Done. Output: $RELATED"
