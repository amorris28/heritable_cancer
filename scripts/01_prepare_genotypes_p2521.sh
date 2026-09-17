#!/usr/bin/env bash
#SBATCH --job-name=prepare_genotypes_p2521
#SBATCH --account=p2521_tsd
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=7792M
#SBATCH --time=04:00:00
#SBATCH --output=prepare_genotypes_p2521.log
#
# One-time, long-running genotype preprocessing for the p2521 (HUSK) cohort:
#   VCF -> BIM -> harmonize PHS601 model -> filter+extract SNPs from VCF -> merge -> score.
# Idempotent: each step is skipped if its output already exists (see scripts/lib/common.sh).
# Submit with: sbatch scripts/01_prepare_genotypes_p2521.sh
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"
load_cohort_config p2521
cd "$REPO_ROOT"

# TSD singularity does not auto-bind paths outside the repo; bind them explicitly.
SINGULARITY_BIND_ARGS=(--bind "$REPO_ROOT" --bind "$GENOTYPE_DIR")

mkdir -p data/interim/vcf_to_bim data/interim/extract_snps data/interim/merge_bfiles data/processed

PHS_FILE=data/external/01_geneticScore.finalModel.snplist.xlsx
MERGED_CSV=data/processed/phs601_${CLUSTER}.csv
MERGED_SNPS=data/processed/phs601_${CLUSTER}.snps
MERGED_BFILE=data/interim/merge_bfiles/${CLUSTER}_phs601
SCORE_FILE=data/processed/${CLUSTER}_phs601.profile

echo "== Step 1: VCF -> BIM (per chromosome, for harmonization only) =="
bim_files=()
for vcf in $VCF_TEMPLATE; do
  base="$(basename "$vcf")"; base="${base%.vcf.gz}"
  out="data/interim/vcf_to_bim/${base}.bim"
  bim_files+=("$out")
  if should_run "$out"; then
    singularity exec "${SINGULARITY_BIND_ARGS[@]}" "$GWAS_IMAGE" plink \
      --vcf "$vcf" \
      --make-just-bim \
      --out "data/interim/vcf_to_bim/${base}"
  fi
done

echo "== Step 2: harmonize PHS601 model per chromosome =="
harmonized_csvs=()
for bim in "${bim_files[@]}"; do
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

echo "== Step 3: merge harmonized per-chr outputs =="
if should_run "$MERGED_CSV" "$MERGED_SNPS"; then
  awk 'FNR==1 && NR!=1{next} {print}' "${harmonized_csvs[@]}" > "$MERGED_CSV"
  awk '{print $6}' "$MERGED_CSV" | tail -n +2 > "$MERGED_SNPS"
fi

echo "== Step 4: filter+extract PHS601 SNPs directly from VCF (per chromosome) =="
filtered_bases=()
for vcf in $VCF_TEMPLATE; do
  base="$(basename "$vcf")"; base="${base%.vcf.gz}"
  out="data/interim/extract_snps/${base}.phs601"
  filtered_bases+=("$out")
  if should_run "${out}.bed" "${out}.bim" "${out}.fam"; then
    singularity exec "${SINGULARITY_BIND_ARGS[@]}" "$GWAS_IMAGE" plink2 \
      --vcf "$vcf" \
      --extract "$MERGED_SNPS" \
      --make-bed \
      --out "$out" \
      --memory 4000
  fi
done

echo "== Step 5: merge filtered bfiles across chromosomes =="
if should_run "${MERGED_BFILE}.bed" "${MERGED_BFILE}.bim" "${MERGED_BFILE}.fam"; then
  merge_list="data/interim/merge_bfiles/merge-list.txt"
  printf '%s\n' "${filtered_bases[@]}" > "$merge_list"
  singularity exec "${SINGULARITY_BIND_ARGS[@]}" "$GWAS_IMAGE" plink \
    --merge-list "$merge_list" \
    --make-bed \
    --out "$MERGED_BFILE"
fi

echo "== Step 6: score merged bfile with PHS601 model =="
if should_run "$SCORE_FILE"; then
  singularity exec "${SINGULARITY_BIND_ARGS[@]}" "$GWAS_IMAGE" plink \
    --bfile "$MERGED_BFILE" \
    --score "$MERGED_CSV" 6 7 5 header sum \
    --out "data/processed/${CLUSTER}_phs601"
fi

echo "Done. Output: $SCORE_FILE"
