#!/usr/bin/env bash
# Dev utility: move a cohort's existing pipeline outputs out of data/public/
# (or, in "manuscript" mode, out of manuscript/figures_and_tables/) into a
# tagged baseline directory, so they survive a rerun of scripts/05-12 and can
# be diffed against freshly regenerated outputs. Use this before re-running
# the pipeline after a code change you want to regression-test (e.g. the
# hazrd -> raw coxph rewrite).
#
# Usage:
#   scripts/backup_outputs_for_diff.sh <p212|p2521> <tag>
#   scripts/backup_outputs_for_diff.sh manuscript <tag>
#
# Example:
#   scripts/backup_outputs_for_diff.sh p212 pre_coxph_refactor
#   # rerun scripts/05-08 for p212 — should_run will see the outputs as
#   # missing (they've been moved) and regenerate them without needing --force
#   scripts/05_fit_models.sh p212
#   scripts/06_plot_figures.sh p212
#   scripts/07_incidence_analysis.sh p212
#   scripts/08_wgs_analysis.sh
#   # then compare:
#   Rscript scripts/compare_outputs.R p212 pre_coxph_refactor
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

MODE="${1:?Usage: $0 <p212|p2521|manuscript> <tag>}"
TAG="${2:?Usage: $0 <p212|p2521|manuscript> <tag>}"
cd "$REPO_ROOT"

if [[ "$MODE" == "manuscript" ]]; then
  # Output of scripts/12_make_figures_and_tables.sh (run locally, not on TSD).
  SRC_DIR="manuscript/figures_and_tables"
  BASELINE_DIR="manuscript/figures_and_tables/_baseline/${TAG}"
  [[ -d "$SRC_DIR" ]] || { echo "No such directory: $SRC_DIR — nothing to back up." >&2; exit 0; }
  mkdir -p "$BASELINE_DIR"
  shopt -s nullglob
  files=("$SRC_DIR"/*)
  shopt -u nullglob
  moved=0
  for f in "${files[@]}"; do
    [[ -f "$f" ]] || continue          # skips directories, including _baseline/
    mv "$f" "$BASELINE_DIR/"
    moved=$((moved + 1))
  done
  echo "Backed up $moved file(s) from $SRC_DIR/ to $BASELINE_DIR/"
  exit 0
fi

COHORT="$MODE"
[[ "$COHORT" == "p212" || "$COHORT" == "p2521" ]] || {
  echo "Usage: $0 <p212|p2521|manuscript> <tag>" >&2
  exit 1
}

# Exact outputs of scripts/05_fit_models.sh, 06_plot_figures.sh, and
# 07_incidence_analysis.sh (see each script's `outputs=(...)` array). Deliberately
# not a glob over data/public/${COHORT}_* — that would also sweep up
# scripts/02_preprocess_metadata.sh outputs (age_summary.tsv etc.), which this
# refactor didn't touch and don't need to be regenerated/compared.
FILES=(
  "${COHORT}_phs_zph_tidy.tsv"
  "${COHORT}_any_pca_performance.tsv"
  "${COHORT}_metastatic_performance.tsv"
  "${COHORT}_pca_death_performance.tsv"
  "${COHORT}_sample_size.tsv"
  "${COHORT}_coxph.tsv"
  "${COHORT}_any_pca_abs_risk.rds"
  "${COHORT}_any_pca_risk_table.rds"
  "${COHORT}_metastatic_abs_risk.rds"
  "${COHORT}_metastatic_risk_table.rds"
  "${COHORT}_zph_tt.rds"
  "${COHORT}_zph_tt_aggr.rds"
  "${COHORT}_cumulative_comparison.rds"
  "${COHORT}_stratified_incidence.rds"
  "${COHORT}_cumulative_comparison_metastatic.rds"
  "${COHORT}_stratified_incidence_metastatic.rds"
)
# scripts/08_wgs_analysis.sh is p212-only.
if [[ "$COHORT" == "p212" ]]; then
  FILES+=(
    p212_wgs_summary.tsv
    p212_fisher.tsv
    p212_wgs_phs_performance_95.tsv
    p212_binomial_phs_hoxb13.tsv
    p212_wgs.tsv
    p212_wgs_cindex.tsv
    p212_wgs_phs_abs_risk.rds
    p212_wgs_phs_risk_table.rds
    p212_wgs_hoxb13_abs_risk.rds
    p212_wgs_hoxb13_risk_table.rds
  )
fi

BASELINE_DIR="data/public/_baseline/${TAG}"
mkdir -p "$BASELINE_DIR"

moved=0
for f in "${FILES[@]}"; do
  if [[ -f "data/public/$f" ]]; then
    mv "data/public/$f" "$BASELINE_DIR/"
    moved=$((moved + 1))
  else
    echo "note: data/public/$f not present, skipping" >&2
  fi
done

echo "Backed up $moved file(s) to ${BASELINE_DIR}/"
echo "Now rerun scripts/05-08 for $COHORT (outputs are absent, so should_run will proceed without --force)."
