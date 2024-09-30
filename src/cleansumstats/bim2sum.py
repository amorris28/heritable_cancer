"""
Creates meta and sumstats files for the cleansumstats pipeline based on the input bim file.
"""

import os
import argparse
import yaml
import pandas as pd

parser = argparse.ArgumentParser(description="Based on the input bim file, creates meta and sumstats files for the cleansumstats pipeline.")
parser.add_argument('--bfile', metavar='meta',help='Absolute file path of the original bed/fam/bim files.', required=True)
parser.add_argument('--out', metavar='out',help='Temporary output folder of this script for the output bim files (.yaml/.gz). Output files (.yaml/.gz) are then used by the cleansumstat pipeline.', required=True)
args = parser.parse_args()

# %%
bfile_dir_in=args.bfile
bfile=os.path.basename(bfile_dir_in)
out_dir = args.out+'/'

# %%  read in bim
print("Reading input bim file '" + bfile_dir_in+'.bim' + "'.")
bim = pd.read_csv(bfile_dir_in+'.bim', sep=None, header=None, names=["CHR","RSID","POSITION","BP","ALT","REF"], engine='python')
bim["BETA"]=1

# %%  create meta file
meta_keys=['cleansumstats_metafile_kind','path_sumStats','stats_Model','stats_neglog10P','col_CHR','col_SNP','col_Direction','col_POS','col_EffectAllele','col_OtherAllele','col_BETA']
meta_items= ['minimal',bfile+'.gz','linear',False] + bim.columns.to_list()
meta=dict(zip(meta_keys,meta_items))

# %% write bim.gz/yaml
print("Writing processed '"+ bfile + ".gz/.yaml' to working directory '" + out_dir +"'.")
bim.to_csv(out_dir + '/' + bfile +'.gz', sep='\t', index=False, compression='gzip')
with open(out_dir + '/' + bfile + '.yaml', 'w') as file:
    yaml.dump(meta, file, sort_keys=True)
print("Done!")