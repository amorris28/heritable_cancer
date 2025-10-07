#!/bin/bash
set -euo pipefail

for chr in {1..22}; do
  echo "Processing chromosome $chr..."
  plink2 --pfile /data/all_hg38 vzs \
	 --rm-dup force-first \
         --chr $chr \
         --make-bed \
         --max-alleles 2 \
         --out /data/chr${chr} \
	 --memory 6000
done

echo "✅ Done!"

