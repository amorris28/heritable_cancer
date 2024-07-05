"""
Removes a user defined region and NaNs. Drops the MHC region by default.
"""

import argparse
import pandas as pd

# parse input arguments
parser = argparse.ArgumentParser(description="Removes a user defined region and NaNs. Drops the MHC region by default.")
parser.add_argument('--input', help='Input path of sumstats.', required=True)
parser.add_argument('--output', help='Output path/name of sumstats.', required=True)
parser.add_argument('--chr', default=6,help='Chromosome number.', required=False)
parser.add_argument('--start', default=26e6, help='Start base pair position.', required=False)
parser.add_argument('--end', default=34e6, help='End base pair position.')
args = parser.parse_args()

chr=str(args.chr)
excl=(int(args.start),int(args.end))

df=pd.read_csv(args.input, sep='\t', dtype=str)
df['POS']=df['POS'].astype(int)
region_df=df[(df['CHR'] == chr) & (df['POS']>=excl[0]) & (df['POS'] <= excl[1])]
print('Excluding {} out of {} total SNPs in chr {} from {} to {}.'.format(region_df.shape[0], df.shape[0], chr, excl[0], excl[1]))
df_exc=df[~df.isin(region_df)].dropna()
df_exc['POS']=df_exc['POS'].astype(int)
df_exc.to_csv(args.output, sep='\t', index=False)
print("Exclusion done!")