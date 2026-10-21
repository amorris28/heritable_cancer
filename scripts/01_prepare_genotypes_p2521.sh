#!/usr/bin/env bash
#SBATCH --job-name=prepare_genotypes_p2521
#SBATCH --account=p2521
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=7792M
#SBATCH --time=04:00:00
#SBATCH --output=prepare_genotypes_p2521.log
#
# One-time, long-running genotype preprocessing for the p2521 (HUSK) cohort:
#   imputed bfile -> harmonize PHS601 model -> extract SNPs -> score.
# Idempotent: each step is skipped if its output already exists (see scripts/lib/common.sh).
# Submit with: sbatch scripts/01_prepare_genotypes_p2521.sh
set -euo pipefail
# sbatch copies this script into a spool dir, so BASH_SOURCE can't locate the repo; use SLURM_SUBMIT_DIR instead.
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  SCRIPT_DIR="${SLURM_SUBMIT_DIR}/scripts"
else
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
source "${SCRIPT_DIR}/lib/common.sh"
load_cohort_config p2521
cd "$REPO_ROOT"

# TSD singularity does not auto-bind paths outside the repo; bind them explicitly.
SINGULARITY_BIND_ARGS=(--bind "$REPO_ROOT" --bind "$GENOTYPE_DIR")

mkdir -p data/interim/merge_bfiles data/processed

PHS_FILE=data/external/01_geneticScore.finalModel.snplist.xlsx
MERGED_CSV=data/processed/phs601_${CLUSTER}.csv
MERGED_SNPS=data/processed/phs601_${CLUSTER}.snps
MERGED_BFILE=data/interim/merge_bfiles/${CLUSTER}_phs601
SCORE_FILE=data/processed/${CLUSTER}_phs601.profile

echo "== Step 1: harmonize PHS601 model against the imputed bim =="
harmonized_csvs=()
for bim in $BIM_TEMPLATE; do
  base="$(basename "$bim" .bim)"
  out="data/processed/phs601_${CLUSTER}_${base}.csv"
  harmonized_csvs+=("$out")
  if should_run "$out"; then
    singularity exec "${SINGULARITY_BIND_ARGS[@]}" "$PYTHON3_IMAGE" python3 src/harmonize_phs.py \
      --phs_file "$PHS_FILE" \
      --bim_file "$bim" \
      --build "$BUILD" \
      --out "$out"
  fi
done

echo "== Step 2: merge harmonized outputs =="
if should_run "$MERGED_CSV" "$MERGED_SNPS"; then
  awk 'FNR==1 && NR!=1{next} {print}' "${harmonized_csvs[@]}" > "$MERGED_CSV"
  awk '{print $6}' "$MERGED_CSV" | tail -n +2 > "$MERGED_SNPS"
fi

echo "== Step 3: extract PHS601 SNPs from the bfile =="
if should_run "${MERGED_BFILE}.bed" "${MERGED_BFILE}.bim" "${MERGED_BFILE}.fam"; then
  for bim in $BIM_TEMPLATE; do
    singularity exec "${SINGULARITY_BIND_ARGS[@]}" "$GWAS_IMAGE" plink2 \
      --bfile "${bim%.bim}" \
      --extract "$MERGED_SNPS" \
      --make-bed \
      --out "$MERGED_BFILE" \
      --memory 4000
  done
fi

echo "== Step 4: score merged bfile with PHS601 model =="
if should_run "$SCORE_FILE"; then
  singularity exec "${SINGULARITY_BIND_ARGS[@]}" "$GWAS_IMAGE" plink \
    --bfile "$MERGED_BFILE" \
    --score "$MERGED_CSV" 6 7 5 header sum \
    --out "data/processed/${CLUSTER}_phs601"
fi

echo "Done. Output: $SCORE_FILE"
