#!/usr/bin/env bash
# Rare-variant (WGS) analysis — p212 only (HOXB13/BRCA2 etc. carrier screening vs PHS).
# Usage: scripts/08_wgs_analysis.sh [--force]
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

[[ "${1:-}" == "--force" ]] && FORCE=1
load_cohort_config p212
cd "$REPO_ROOT"

if [[ "$HAS_WGS" != "1" ]]; then
  echo "WGS analysis is only defined for p212; nothing to do." >&2
  exit 0
fi

PHENOTYPES=data/processed/phenotypes.tsv
[[ -f "$PHENOTYPES" ]] || { echo "Missing $PHENOTYPES — run scripts/02_preprocess_metadata.sh p212 first" >&2; exit 1; }

outputs=(
  data/public/p212_wgs_summary.tsv
  data/public/p212_fisher.tsv
  data/public/p212_wgs_phs_performance_95.tsv
  data/public/p212_binomial_phs_hoxb13.tsv
  data/public/p212_wgs.tsv
  data/public/p212_wgs_cindex.tsv
  data/public/p212_wgs_cuminc_at_age.tsv
  data/public/p212_wgs_phs_abs_risk.rds
  data/public/p212_wgs_phs_risk_table.rds
  data/public/p212_wgs_hoxb13_abs_risk.rds
  data/public/p212_wgs_hoxb13_risk_table.rds
)
should_run "${outputs[@]}" || exit 0

run_rscript src/fit-models-plot-figures-wgs.R "$PHENOTYPES" TRUE p212 "$BOOTSTRAP_ITERS"
require_files p212_wgs_summary.tsv p212_fisher.tsv p212_wgs_phs_performance_95.tsv \
  p212_binomial_phs_hoxb13.tsv p212_wgs.tsv p212_wgs_cindex.tsv p212_wgs_cuminc_at_age.tsv \
  p212_wgs_phs_abs_risk.rds p212_wgs_phs_risk_table.rds \
  p212_wgs_hoxb13_abs_risk.rds p212_wgs_hoxb13_risk_table.rds

mv p212_wgs_summary.tsv p212_fisher.tsv p212_wgs_phs_performance_95.tsv \
   p212_binomial_phs_hoxb13.tsv p212_wgs.tsv p212_wgs_cindex.tsv p212_wgs_cuminc_at_age.tsv \
   p212_wgs_phs_abs_risk.rds p212_wgs_phs_risk_table.rds \
   p212_wgs_hoxb13_abs_risk.rds p212_wgs_hoxb13_risk_table.rds \
   data/public/

echo "Done. Outputs: ${outputs[*]}"
