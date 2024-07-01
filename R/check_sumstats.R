library(tidyverse)

neff = function(ncases, ncontrols) { 4 / (1 / ncases + 1 / ncontrols)}

sumstats = read_tsv('PGS000067_hmPOS_GRCh37.txt', comment = "#")

head(sumstats)

# number of SNPs in PRS
nrow(sumstats)

# rename columns and drop unnecessary columns
sumstats %>% 
  select(CHR = chr_name, 
         SNP = rsID,
         BP = chr_position,
         A1 = effect_allele,
         A2 = other_allele,
         OR = effect_weight
         ) ->
sumstats_cleaned 

# ambiguous SNPs A/T C/G
sumstats_cleaned %>% 
  filter((A1 %in% c('C', 'G') & A2 %in% c('C', 'G')) | 
           (A1 %in% c('A', 'T') & A2 %in% c('A', 'T')))

# remove C -> G, G -> C, A -> T, and T-> A
sumstats_cleaned %>% 
  filter(!((A1 %in% c('C', 'G') & A2 %in% c('C', 'G')) | 
           (A1 %in% c('A', 'T') & A2 %in% c('A', 'T')))) ->
  sumstats_cleaned
nrow(sumstats_cleaned)

# non-ACTG SNPs
sumstats_cleaned %>% 
  filter(!(A1 %in% c('A', 'C', 'T', 'G')) | !(A2 %in% c('A', 'C', 'T', 'G')))

# SNPs without an rs#
sumstats_cleaned %>% 
  filter(!str_starts(SNP, 'rs') | is.na(SNP))

write_tsv(sumstats_cleaned, 'PGS000067_hmPOS_GRCh37_cleaned.txt')
