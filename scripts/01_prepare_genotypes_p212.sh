#!/usr/bin/env bash
#SBATCH --job-name=prepare_genotypes_p212
#SBATCH --account=p212_tsd
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=7792M
#SBATCH --time=04:00:00
#SBATCH --output=prepare_genotypes_p212.log
#
# One-time, long-running genotype preprocessing for the p212 (NFCS/MoBA) cohort:
#   BIM files (already imputed) -> harmonize PHS601 model -> filter SNPs -> merge -> score -> HOXB13.
# Idempotent: each step is skipped if its output already exists (see scripts/lib/common.sh).
# Submit with: sbatch scripts/01_prepare_genotypes_p212.sh
# Force recomputation of a step by removing its output file(s), or FORCE=1 sbatch ...
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"
load_cohort_config p212
cd "$REPO_ROOT"

# TSD singularity does not auto-bind paths outside the repo; bind them explicitly.
SINGULARITY_BIND_ARGS=(--bind "$REPO_ROOT" --bind "$GENOTYPE_DIR")

mkdir -p data/interim/filter_bfiles data/interim/merge_bfiles data/processed

PHS_FILE=data/external/01_geneticScore.finalModel.snplist.xlsx
MERGED_CSV=data/processed/phs601_${CLUSTER}.csv
MERGED_SNPS=data/processed/phs601_${CLUSTER}.snps
MERGED_BFILE=data/interim/merge_bfiles/${CLUSTER}_phs601
SCORE_FILE=data/processed/${CLUSTER}_phs601.profile
HOXB13_FILE=data/processed/tmp_rs138213197.raw

echo "== Step 1: harmonize PHS601 model per chromosome =="
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

echo "== Step 2: merge harmonized per-chr outputs =="
if should_run "$MERGED_CSV" "$MERGED_SNPS"; then
  awk 'FNR==1 && NR!=1{next} {print}' "${harmonized_csvs[@]}" > "$MERGED_CSV"
  awk '{print $6}' "$MERGED_CSV" | tail -n +2 > "$MERGED_SNPS"
fi

echo "== Step 3: filter bfiles to PHS601 SNPs (per chromosome) =="
filtered_bases=()
for bim in $BIM_TEMPLATE; do
  base="$(basename "$bim" .bim)"
  dir="$(dirname "$bim")"
  out="data/interim/filter_bfiles/${base}.phs601"
  filtered_bases+=("$out")
  if should_run "${out}.bed" "${out}.bim" "${out}.fam"; then
    singularity exec "${SINGULARITY_BIND_ARGS[@]}" "$GWAS_IMAGE" plink2 \
      --bfile "${dir}/${base}" \
      ${PSAM_FLAG} \
      --extract "$MERGED_SNPS" \
      --make-bed \
      --out "$out"
  fi
done

echo "== Step 4: merge filtered bfiles across chromosomes =="
if should_run "${MERGED_BFILE}.bed" "${MERGED_BFILE}.bim" "${MERGED_BFILE}.fam"; then
  merge_list="data/interim/merge_bfiles/merge-list.txt"
  printf '%s\n' "${filtered_bases[@]}" > "$merge_list"
  singularity exec "${SINGULARITY_BIND_ARGS[@]}" "$GWAS_IMAGE" plink \
    --merge-list "$merge_list" \
    --make-bed \
    --out "$MERGED_BFILE"
fi

echo "== Step 5: score merged bfile with PHS601 model =="
if should_run "$SCORE_FILE"; then
  singularity exec "${SINGULARITY_BIND_ARGS[@]}" "$GWAS_IMAGE" plink \
    --bfile "$MERGED_BFILE" \
    --score "$MERGED_CSV" 6 7 5 header sum \
    --out "data/processed/${CLUSTER}_phs601"
fi

echo "== Step 6: extract HOXB13 (rs138213197) genotype =="
if should_run "$HOXB13_FILE"; then
  singularity exec "${SINGULARITY_BIND_ARGS[@]}" "$GWAS_IMAGE" plink \
    --bfile "$MERGED_BFILE" \
    --snp rs138213197 \
    --recode A \
    --out data/processed/tmp_rs138213197
fi

echo "Done. Outputs: $SCORE_FILE, $HOXB13_FILE"
