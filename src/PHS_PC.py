#!/usr/bin/env python3
import argparse
import pandas as pd
from pathlib import Path

def merge_by_chr_bp_and_flip_alleles(bim, phs, chrcol, bpcol):
    # NB - do not consider strand flips, and retain ambiguous SNPs
    # NB! There are duplicates on CHR:BP, but it's ok as long as we don't have duplicates after merging with phs file

    print(f'{len(phs)} SNPs in original model')

    #assert 0==bim['SNP'].duplicated().sum()
    if 'GP' in bim: del bim['GP']

    phs = pd.merge(phs[[chrcol, bpcol, 'Effect', 'Ref', 'Beta']].rename(columns={chrcol:'CHR',bpcol:'BP'}), bim, on=['CHR','BP'], how='left')
    phs=phs[~phs['SNP'].isnull()]
    print(f'{len(phs)} SNPs after merging on CHR:BP')
    
    idx = (phs['A1']==phs['Effect']) & (phs['A2']==phs['Ref'])
    idx2 = (phs['A1']==phs['Ref']) & (phs['A2']==phs['Effect']) 
    
    phs.loc[idx2, 'Beta'] = -phs.loc[idx2, 'Beta'] 
    phs.loc[idx2, ['Effect', 'Ref']] = phs.loc[idx2, ['Ref', 'Effect']].values
    print(f'Flip A1<->A2 and invert sign of BETA for {idx2.sum()} SNPs')
    
    phs = phs[idx | idx2]
    assert (phs['Effect'] == phs['A1']).all()
    assert (phs['Ref'] == phs['A2']).all()
    print(f'{len(phs)} SNPs after matching alleles')
    
    assert not bim[bim['SNP'].isin(phs['SNP'].values)].duplicated().any()  # check no duplicated SNPs
    #assert not phs['SNP'].duplicated().any()  # check no duplicated SNPs
    return phs

def main():
    parser = argparse.ArgumentParser(description="Harmonize model sumstats with genotype .bim files.")
    parser.add_argument("--phs_file", required=True, help="Path to model sumstats Excel file.")
    parser.add_argument("--bim_file", required=True, help="BIM file")
    parser.add_argument("--build", choices=["37","38"], default="37", help="Genome build for position column.")
    parser.add_argument("--out", required=True, help="Output file path (.csv).")
    args = parser.parse_args()

    bpcol = "Position (GRCh37)" if args.build == "37" else "Position (GRCh38)"
    chrcol = "Chromosome"

    print(f"Loading model file: {args.phs_file}")
    phs = pd.read_excel(args.phs_file)

    bim_path = Path(args.bim_file)
    if not bim_path.exists():
        raise FileNotFoundError(f"⚠️ BIM file not found: {bim_path}")

    bim = pd.read_csv(bim_path, sep='\t', header=None, 
                      names=['CHR','SNP','GP','BP','A1','A2'])
    print(f"{len(bim)} SNPs loaded from {bim_path}")

    phs_harmonized = merge_by_chr_bp_and_flip_alleles(bim, phs, chrcol, bpcol)
    phs_harmonized.to_csv(args.out, sep='\t', index=False)
    print(f"✅ Wrote harmonized file to {args.out}")

if __name__ == "__main__":
    main()

