#!/usr/bin/env bash
set -euo pipefail
mkdir -p data/interim/king
../bin/king -b ../HUSK_v240619/husk_all.autosomes.v240619.genomatch.bed \
  --related --degree 3 --prefix data/interim/king/HUSK_related_deg3