#!/usr/bin/env python3
"""Convert a KING text lower-triangle kinship matrix to plink2 --make-king
'triangle bin' format (float64, lower triangle, no diagonal).

usage: king_text_to_bin.py kinship_results25.king kinship_results25.king.id out_prefix
writes out_prefix.king.bin and out_prefix.king.id
"""
import shutil
import sys

import numpy as np

text_path, id_path, out_prefix = sys.argv[1:4]

with open(id_path) as f:
    ids = [line for line in f if line.strip()]
if ids[0].startswith("#"):
    ids = ids[1:]
n = len(ids)
expected = n * (n - 1) // 2
print(f"samples in id file: {n}; expected values: {expected}")

written = 0
with open(text_path) as src, open(f"{out_prefix}.king.bin", "wb") as out:
    for line in src:
        values = np.array(line.split(), dtype=np.float64)
        if values.size:
            values.tofile(out)
            written += values.size

print(f"values written: {written}")
if written != expected:
    sys.exit(
        f"ERROR: {written} values but {n} ids imply {expected}. "
        "Check that the .king and .king.id files match."
    )

shutil.copyfile(id_path, f"{out_prefix}.king.id")
print(f"wrote {out_prefix}.king.bin ({written * 8} bytes) and {out_prefix}.king.id")
