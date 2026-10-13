#!/usr/bin/env bash
#SBATCH --job-name=cohort_pca
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=01:00:00
#SBATCH --output=cohort_pca_%j.log
set -euo pipefail

# sbatch copies this script into a spool dir, so BASH_SOURCE can't locate the repo.
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

OUT="data/interim/pca"
INPUT_BFILE="$BED_FILE"
PRUNE_PREFIX="$OUT/$COHORT"
PCA_BFILE="$OUT/$COHORT.pruned"

case "$BUILD" in
  37)
    HIGH_LD_FILE="${SCRIPT_DIR}/resources/high-LD-regions-hg19-GRCh37.txt"
    ;;
  38)
    HIGH_LD_FILE="${SCRIPT_DIR}/resources/high-LD-regions-hg38-GRCh38.txt"
    ;;
  *)
    echo "Unsupported genome build for PCA: $BUILD" >&2
    exit 1
    ;;
esac
[[ -f "$HIGH_LD_FILE" ]] || {
  echo "Missing build $BUILD high-LD regions file: $HIGH_LD_FILE" >&2
  exit 1
}

mkdir -p "$OUT"
PLINK="../bin/plink"
FLASHPCA="$(cd ../bin && pwd)/flashpca_x86-64"
[[ -x "$FLASHPCA" ]] || { echo "Missing or non-executable $FLASHPCA" >&2; exit 1; }
PCA_CPUS="${SLURM_CPUS_PER_TASK:-$(nproc 2>/dev/null || echo 1)}"

[[ -n "$INPUT_BFILE" ]] || { echo "BED_FILE is not set in config/$COHORT.env" >&2; exit 1; }
require_files "$INPUT_BFILE.bed" "$INPUT_BFILE.bim" "$INPUT_BFILE.fam"

if should_run "$PRUNE_PREFIX.prune.in"; then
  "$PLINK" \
    --bfile "$INPUT_BFILE" \
    --chr 1-22 --geno 0.05 \
    --indep-pairwise 1000 50 0.05 \
    --exclude range "$HIGH_LD_FILE" \
    --threads "$PCA_CPUS" \
    --out "$PRUNE_PREFIX"
fi

if should_run "$PCA_BFILE.bed" "$PCA_BFILE.bim" "$PCA_BFILE.fam"; then
  "$PLINK" \
    --bfile "$INPUT_BFILE" \
    --extract "$PRUNE_PREFIX.prune.in" \
    --make-bed \
    --threads "$PCA_CPUS" \
    --out "$PCA_BFILE"
fi
require_files "$PCA_BFILE.bed" "$PCA_BFILE.bim" "$PCA_BFILE.fam"

PCA_OUTPUTS=(
  "$OUT/eigenvectors_${CLUSTER}.txt"
  "$OUT/pcs_${CLUSTER}.txt"
  "$OUT/eigenvalues_${CLUSTER}.txt"
  "$OUT/pve_${CLUSTER}.txt"
  "$OUT/loadings_${CLUSTER}.txt"
  "$OUT/meansd_${CLUSTER}.txt"
)
if should_run "${PCA_OUTPUTS[@]}"; then
  (
    cd "$OUT"
    "$FLASHPCA" \
      --bfile "$COHORT.pruned" \
      --outload "loadings_${CLUSTER}.txt" \
      --outmeansd "meansd_${CLUSTER}.txt" \
      --numthreads "$PCA_CPUS" \
      --suffix "_${CLUSTER}.txt"
  )
fi
require_files "${PCA_OUTPUTS[@]}"

echo "Done. Pruned genotype set: $PCA_BFILE"
echo "FlashPCA PCs: $OUT/pcs_${CLUSTER}.txt"
