"""
Performs some preprocessing steps on the meta and sumstat files before executing the cleansumstats pipeline.
Output path needs to be different form input folder to avoid overwriting sumstats and meta file.
"""

import os
import argparse
from argparse import RawTextHelpFormatter
import yaml
import pandas as pd
import numpy as np
import re 

# parse input arguments
parser = argparse.ArgumentParser(description="Preprocessing steps for the cleansumstats pipeline. Make sure that either \n - col_N \n - CaseN and ControlN \n - stats_TotalN \n - stats_caseN and stats_ControlN \n are defined in the metafile ('.yaml'). Labels containinng special characters must be placed in '' (e.g., '#chrom') and will be removed from the processed meta and sumstat file.",formatter_class=RawTextHelpFormatter)
parser.add_argument('--meta', metavar='meta',help='Absolute file path of the meta file (.yaml).', required=True)
parser.add_argument('--out', metavar='out',help='Temporary output folder of this script for the pocessed files (.yaml/sumstats). Needs to be different from the metafile location. Output files (.yaml/sumstats) are then used by the cleansumstat pipeline.', required=True)
parser.add_argument('--sumstat', metavar='sumstat',help='Optional: sumstat file path. Else it is assumed to reside next to the metafile (.yaml).')
parser.add_argument("--sep", default=None, type=str, help="Delimiter choice ('\\t','whitespace', ' ', ',', ';') for reading in sumstats with in pandas.read_csv(). Default: 'sep=None'")
args = parser.parse_args()

# assign input arguments
with open(args.meta, 'r') as file:
    meta = yaml.safe_load(file)
metadir = os.path.dirname(args.meta) # gives the absolute path of the meta directory
sumstatdir = args.sumstat if args.sumstat else metadir + '/' + meta["path_sumStats"]
trait,_ = os.path.splitext(args.meta)
trait = os.path.basename(trait) # trait name = file name of the meta file (without folder)
out = args.out
assert os.path.normpath(metadir) != os.path.normpath(out), 'Output folder must be different form input folder to avoid overwriting the input files (sumstats and meta file)'
sep ='\s+' if args.sep=='whitespace' else args.sep
engine='python' if sep==None else None

# read in sumstats
print("Reading input sumstats '" + sumstatdir + "'.")
sumstat = pd.read_csv(sumstatdir, sep=sep, engine=engine)

# check for cleansumstat_metafile_kind and change "stats_Model" when it is missing in the raw sumstats (needed for producing z col)
meta.update({"cleansumstats_metafile_kind": "minimal"}) if 'cleansumstats_metafile_kind' not in meta else None
meta["stats_Model"]="linear" if 'stats_Model' in meta and meta["stats_Model"]=="linear mixed-model" else meta["stats_Model"]
meta["stats_Model"]="logistic" if 'stats_Model' in meta and meta["stats_Model"]=="logistic mixed-model" else meta["stats_Model"]

# raise an error if a column is defined in meta, but it is missing from the sumstat file this would indicate an invalid meta file; a change in meta file is required
cols = ['col_N', 'col_CaseN', 'col_ControlN'] #,'col_SE','col_POS']
# for col in cols:
    # if (col in meta) and (meta[col] not in sumstat.columns):
#         # raise(ValueError(f'column {meta[col]} is defined by meta[{col}], but is missing from the sumstats file'))
assert np.isin(np.array(list(meta.values()))[np.isin(list(meta.keys()),cols)], sumstat.columns.to_list()).all(), 'The value "col_N", "col_CaseN" or "col_ControlN" is defined in the metafile but is missing/misspelled in the sumstats file!'

# check for N, CaseN and ControlN in sumstats and meta file and update files with proper details:
#if 'col_N' in meta:
#    print("N is in sumstat!")
#elif ('col_CaseN' in meta) and ('col_ControlN' in meta):
#    print("CaseN and ControlN are in sumstat!")
#    meta.update({"col_N": "TotalSample"})
#    sumstat['TotalSample']=np.divide(4, np.divide(1.0, pd.to_numeric(sumstat[meta["col_CaseN"]]).values) + np.divide(1.0, pd.to_numeric(sumstat[meta["col_ControlN"]]).values)) 
#elif 'stats_TotalN' in meta:
#    print("stats_TotalN is in meta file!")
#    meta.update({"col_N": 'TotalSample'})
#    sumstat['TotalSample'] = meta["stats_TotalN"]
#elif ('stats_CaseN' in meta) and ('stats_ControlN' in meta):
#    print("stats_caseN and stats_ControlN are in meta file!")
#    meta.update({"col_CaseN": 'CaseN',"col_ControlN": 'ControlN',"col_N": 'TotalSample'})
#    sumstat['CaseN'], sumstat['ControlN'], sumstat['TotalSample'] = [meta["stats_CaseN"], meta["stats_ControlN"], 4/((1/meta["stats_CaseN"])+(1/meta["stats_ControlN"]))]
#else:
#    raise ValueError("N, CaseN or ControlN are missing in the sumstat and meta files. Please update your files with the corresponding information!")

# clean sumstats header and meta file items accordingly
sumstat=sumstat.rename(index=str,columns=dict(zip(sumstat.columns,[re.sub('[^a-zA-Z0-9 ]', '', i) for i in sumstat.columns])))
sumstat.columns = sumstat.columns.str.lower()  # lower case all column names to reduce typos by user
assert not sumstat.columns.has_duplicates, 'After transforming column names to lower case, there are duplicats!'
sumstat.columns = sumstat.columns.str.strip() # manage white spaces in the header labels
for key, value in meta.items():
    meta.update({key: re.sub(r'[^a-zA-Z0-9 ]', '', str(value).lower().strip()) if type(value) is str else value})
meta.update({"path_sumStats": trait +'.gz'})

# save updated sumstats and meta files
print(f"Writing processed '{trait}.gz/.yaml' to temporary output folder {out}.")
sumstat.to_csv(f"{out}/{trait}.gz", sep='\t', index=False, compression='gzip')
with open(f"{out}/{trait}.yaml", 'w') as file:
    yaml.dump(meta, file, sort_keys=True)
print("Done!")
