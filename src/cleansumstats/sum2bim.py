"""
Reshapes the output bim file of the cleansumstats pipeline back to the original PLINK structure.
"""

import os
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Reshapes the output bim file of the cleansumstats pipeline back to the original PLINK structure.")
parser.add_argument('--bfile', metavar='meta',help='Absolute file path of the original bed/fam/bim files.', required=True)
parser.add_argument('--bim', metavar='meta',help='Absolute file path of the bim file processed with the cleansumstats pipeline.', required=True)
parser.add_argument('--out', metavar='out',help='Output folder.', required=True)
args = parser.parse_args()

# %%
bfile_dir_in=args.bfile
bim_clsum_dir=args.bim
out_dir=args.out+'/'
bfile_out=os.path.basename(bfile_dir_in)+"_cs"

# %% read files
print(f"Reading original/cleansumstat-processed bim file\n {bfile_dir_in}.bim\n {bim_clsum_dir}.bim.")
bim_org = pd.read_csv(bfile_dir_in+'.bim', sep=None, header=None, names=["CHR","RSID","POSITION","BP","ALT","REF"], engine='python')
bim_clsum = pd.read_csv(bim_clsum_dir, sep=None, engine='python')
bim_clsum=bim_clsum.rename(columns={"Direction": "POSITION", "POS": "BP", "EffectAllele": "ALT", "OtherAllele": "REF"})

# %%  Reshape output bim to original PLINK dimension
assert not bim_clsum["0"].duplicated().any(), 'Duplicate indices!'
assert bim_clsum["B"].shape[0] == (bim_clsum["B"]==1).sum() + (bim_clsum["B"]==-1).sum(), 'n(B==1)+n(B==-1)!=n(B)'
print(f"Multi-allelic RSIDs: {', '.join(bim_clsum.loc[bim_clsum['RSID'].duplicated(),'RSID'].to_list())}") if bim_clsum["RSID"].duplicated().any() else None

# check if SNPs have been remapped to another chr.
if bim_clsum['CHR'].value_counts().index.shape[0]>1 and bim_org['CHR'].value_counts().index.shape[0]<=1:
    print('Non-unique Chr identifiers in bim after cleansumstats processing. Some SNPs likely have been remapped to another chr. Respective SNPs will be removed:'); print(bim_clsum['CHR'].value_counts()) 
    chr=bim_clsum['CHR'].value_counts().index[0]
    bim_clsum=bim_clsum.loc[bim_clsum['CHR']==chr]
    bim_clsum=bim_clsum.astype({'CHR': 'int'})
    assert pd.unique(bim_clsum['CHR']).shape[0]==1

# create bim_out with data from bim_clsum and dimension of bim_org
bim_out=bim_org.copy()
bim_out["RSID"]='.'
betap=[bim_clsum["B"]==1, bim_clsum["0"].loc[bim_clsum["B"]==1]-1]
betan=[bim_clsum["B"]==-1, bim_clsum["0"].loc[bim_clsum["B"]==-1]-1]
for col in bim_out.columns:
    if col not in ["REF","ALT"]:
        bim_out.loc[bim_clsum["0"]-1,col]=bim_clsum[col].values
    else:
        bim_out.loc[betap[1],col]=bim_clsum[col].loc[betap[0]].values # beta==1: overwrite REF/ALT with data from cleansumstats
        if col=="REF": # beta==-1
            bim_out.loc[betan[1],col]=bim_clsum["ALT"].loc[betan[0]].values # place ALT/A1 from cleansumstats in REF/A2 in bim_out to flip bed with plink
        else:
            bim_out.loc[betan[1],col]=bim_clsum["REF"].loc[betan[0]].values # place REF/A2 from cleansumstats in ALT/A1 in bim_out to flip bed with plink   

bim_out.loc[bim_out["RSID"]!='.',"RSID"]= bim_out.loc[bim_out["RSID"]!='.',"RSID"] + '_' + pd.Series(bim_out.loc[bim_out["RSID"]!='.'].index).astype(str).to_list() # add line suffix from RSIDS for multi-allelic SNPs 
bim_out.loc[bim_out["RSID"]=='.','BP']=-bim_out.loc[bim_out["RSID"]=='.','BP'] # will be ignored/removed by plink

# for b==-1 check if SNPs were rotated correctly 
assert (bim_out.loc[betan[1],["ALT","REF"]].values==bim_clsum[["REF","ALT"]].loc[betan[0]].values).all(), 'bim_out[ALT,REF]!=bim_clsum[REF,ALT] for beta==-1!'

# %% export files
print(f"Writing processed '{bim_clsum_dir}' to working directory '{out_dir}'.")
bim_out.to_csv(f"{out_dir}{bfile_out}_tmp.bim",sep="\t",header=None,index=False)
bim_out.loc[betan[1],["RSID","REF"]].to_csv(f"{out_dir}forceALT_betan_{bfile_out}.txt",sep="\t",header=None,index=False)
print("Done!")