import pysam
import pandas as pd

def merge_by_chr_bp_and_flip_alleles(bim, phs, chrcol, bpcol):
    # NB - do not consider strand flips, and retain ambiguous SNPs
    # NB! There are duplicates on CHR:BP, but it's ok as long as we don't have duplicates after merging with phs file

    print(f'{len(phs)} SNPs in original model')

    assert 0==bim['SNP'].duplicated().sum()
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
    
    assert not phs['SNP'].duplicated().any()  # check no duplicated SNPs
    return phs

phs601=pd.read_excel('data/external/01_geneticScore.finalModel.snplist.xlsx')
# bpcol = 'Position (GRCh37)'   # for p212 familiar cancer
bpcol = 'Position (GRCh38)'  # for HUSK

# phs601_chr = []
# for chri in range(1, 23):
#     bim_fname=f'/ess/p212/cluster/moba_qc_imputation/BA/release6/chr{chri}.imputed.bim'
#     bim=pd.read_csv(bim_fname, sep='\t', header=None, names=['CHR','SNP','GP','BP','A1','A2'])
#     print(f'{len(bim)} SNPs loaded from {bim_fname} file')

#     phs601_chr.append(merge_by_chr_bp_and_flip_alleles(bim, phs601, 'Chromosome', bpcol))
#     del bim



phs601_chr = []
for chri in range(1, 23):
    vcf_fname=f'/ess/p2521/data/durable/HuskMent/genotypes/imputed/v240619/vcf/husk_all.imputed.chr{chri}.vcf.gz'
    records = []
    with pysam.VariantFile(vcf_fname) as vcf:
        for i, rec in enumerate(vcf):
            records.append([rec.chrom, rec.pos, rec.id, rec.ref, ','.join(rec.alts)])
            if i % 1_000_000 == 0 and i > 0:
                print(f"{i:,} records processed...")  # optional progress
            # optionally: stop early for testing
            # if i > 9: break

    vcf = pd.DataFrame(records, columns=["CHR","BP","SNP","A1","A2"])
    vcf["CHR"] = vcf["CHR"].astype(int)
    print(f'{len(vcf)} SNPs loaded from {vcf_fname} file')

    phs601_chr.append(merge_by_chr_bp_and_flip_alleles(vcf, phs601, 'Chromosome', bpcol))
    del vcf

phs601 = pd.concat(phs601_chr)
# phs601.to_csv('phs601_heritable_cancer.csv', sep='\t', index=False)
phs601.to_csv('data/interim/phs601_husk.csv', sep='\t', index=False)
