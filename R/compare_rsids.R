library(tidyverse)

directory = "derived_data/"
files = c("rsids_seibert", "rsids_huynle2021", "rsids_huynle2022", "rsids_karunamuni")

sumstat_rsids = lapply(files, function(x) read_tsv(paste0(directory, x), col_names = "RSIDs"))

names(sumstat_rsids) = files

rsids = read_lines("derived_data/RSIDS_ANTE")

sumstat_rsids2 = lapply(sumstat_rsids, function(x) mutate(x, in_antegenes = RSIDs %in% rsids) %>% 
                          filter(!in_antegenes))
bind_rows(sumstat_rsids2, .id = "model") %>% 
  write_tsv("derived_data/missing_rsids_in_genotype_data")

sumstats = read_tsv('derived_data/sumstats/PRACTICAL_PCa_Seibert2018.GRCh37.gz')

sumstats[sumstats$RSID %in% sumstat_rsids2$rsids_seibert$RSIDs, ]
