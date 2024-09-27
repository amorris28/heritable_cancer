#!/usr/bin/env python3

import pandas as pd
import re

colnames = ["FID", "IID", "FA", "MO", "SEX", "PHE"]
fam = pd.read_csv('/cluster/projects/p212/antenor_pipeline/antegenes_imputed_allchr-temporary.fam', names = colnames, sep = '\t', header = None)

fam = fam.iloc[:, 0:2]
new_fam = fam.copy()
def fix_ids(x):
    return re.sub(r"(.*?_.*?)_.*", "\\1", x)

new_fam['IID'] = new_fam["IID"].apply(fix_ids)

output = pd.concat([fam, new_fam], axis = 1)
output.to_csv('data/processed/update-antegenes-ids.tsv', sep = "\t", index = False, header = False, na_rep = "NA")
