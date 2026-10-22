#!/usr/bin/env bash
# Internal absolute-risk calibration of PHS at screening ages (quintiles).
# Usage: scripts/12_calibration.sh <p212|p2521> [--force]
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

COHORT="${1:?Usage: $0 <p212|p2521> [--force]}"
[[ "${2:-}" == "--force" ]] && FORCE=1
load_cohort_config "$COHORT"
cd "$REPO_ROOT"

PHENOTYPES=data/processed/phenotypes.tsv
[[ -f "$PHENOTYPES" ]] || { echo "Missing $PHENOTYPES — run scripts/02_preprocess_metadata.sh $CLUSTER first" >&2; exit 1; }

OUTPUT="data/public/${CLUSTER}_calibration_by_age.tsv"
should_run "$OUTPUT" || exit 0

run_rscript src/calibration.R "$PHENOTYPES" "$CLUSTER"
for f in calibration_by_age.tsv calibration_slope_intercept.tsv calibration.png; do
  require_files "${CLUSTER}_${f}"
  mv "${CLUSTER}_${f}" "data/public/${CLUSTER}_${f}"
done
[[ -f "${CLUSTER}_calibration_delayed_entry.png" ]] && mv "${CLUSTER}_calibration_delayed_entry.png" "data/public/"

echo "Done. Output: $OUTPUT"
