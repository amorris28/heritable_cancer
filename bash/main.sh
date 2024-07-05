
# Pre-process steps: 
# 1. Download sumstats
# 2. Create yaml files (Ex: https://github.com/amorris28/GWAS_SUMSTAT/blob/master/cleansumstats_meta/PGC_SCZ_0518_EUR/PGC_SCZ_0518_EUR.yaml)
# 3. Remove leading lines with comments (#)
# 4. Import sumstats + yaml + BIM file(s) to p33

gzcat raw_data/sumstats/PGS000067_hmPOS_GRCh37.txt.gz \
  | sed '/^#/d' \
  | gzip \
  > derived_data/sumstats/PRACTICAL_PCa_Seibert2018.txt.gz

########################################
# On p33

# Run cleansumstats using JF script from the GWAS_SUMSTATS github
# Customizations:
# Comment out the N section of prep_sum.py
# 
# Can call cleansumstats directly, but I find the wrapper makes it easier
# to work with it on TSD 
# https://github.com/amorris28/GWAS_SUMSTAT/tree/custom/JF
sbatch run_cleansumstats_bulk.job
sbatch run_cleansumstats_bim_bulk.job

# Export {${BFILE_OUT}_tmp.bim,forceALT_betan_${BFILE_OUT}.txt from p33
# Export sumstats.GRCh37.gz
# Import to the local p project for your B files (e.g., p212)

#########################################
# On p212

# Identify cases and controls and sample IDs
module load R/4.2.1-foss-2022a
Rscript identify_cases_controls.R


# For the antegenes/1KG imputation Bfiles:
# Update the Sample IDs since the fid and iid's have been created incorrectly
Rscript update_ids_antegenes.R

###########################################
# On local PC (or p212)
# Should I first identify SNPs in LD with missing SNPs from GWAS?
# 1. Check for RSIDs from sumstats in imputed (rotated) genotype file
# 2. Can I look for LD based on Chr and location even if the RSID is not in BIM?
# 3. Check dbsnp or something for RSIDs in LD
# RSIDs in cleaned sumstats file
SUMSTATS=derived_data/sumstats/
gzcat ${SUMSTATS}PRACTICAL_PCa_Seibert2018.GRCh37.gz | awk 'NR > 1 { print $3 }' > derived_data/seibert_rsids
gzcat ${SUMSTATS}PRACTICAL_PCa_Huynh-Le2021.GRCh37.gz | awk 'NR > 1 { print $3 }' > derived_data/huynle2021_rsids
gzcat ${SUMSTATS}PRACTICAL_PCa_Huynh-Le2022.GRCh37.gz | awk 'NR > 1 { print $3 }' > derived_data/huynle2022_rsids
gzcat ${SUMSTATS}PRACTICAL_PCa_Karunamuni2021.GRCh37.gz | awk 'NR > 1 { print $3 }' > derived_data/karunamuni_rsids

# RSIDs in genotype file (BIM)
#awk '{ print $2 }' /cluster/projects/p212/moba_qc_imputation/BA/release6/moba_imputed_all.bim > MOBA_RSIDS
awk '{ print $2 }' moba_imputed_all_cs.bim > MOBA_RSIDS
awk '{ print $2 }' antegenes_imputed_allchr-temporary_cs.bim > ANTE_RSIDS

# Look up these RSIDs on LDproxy including all European populations
# rs10905266 instead of rs201305416
# rs1010 is not biallelic
# rs1351836062 is not in reference panel

# Check for individual RSIDs
awk -F'\t' '$2 == "rs201305416" {print}' /tsd/p212/home/p212-andrewmo/out/antegenes_imputed_allchr-temporary_cs.bim

# Use matrix multiplication to calcule PHS values for each combinationg of
# model and genotype
for VARIABLE in PRACTICAL_PCa_Seibert2018.GRCh37.gz PRACTICAL_PCa_Huynh-Le2021.GRCh37.gz PRACTICAL_PCa_Huynh-Le2022.GRCh37.gz PRACTICAL_PCa_Karunamuni2021.GRCh37.gz 
do
  zcat $VARIABLE | awk '{ print $3 }' | tail -n +2 >> phs_rsids_test
done
echo "rs10905266" >> phs_rsids_test

# Rotate Alleles


# Calculate PHS scores
# Should add a step before this script to swap in the SNPs in LD
Rscript R/calculate_phs.R sumstats_list.txt out/antegenes_imputed_allchr-temporary_cs

# Fit coxph and plot K-M curves
# Run RK_get_perf.R to get HRs
sbatch bash/phs.batch

