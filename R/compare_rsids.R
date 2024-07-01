library(tidyverse)

directory = "derived_data/"
files = c("seibert_rsids", "huynle2021_rsids", "huynle2022_rsids", "karunamuni_rsids")

sumstat_rsids = lapply(files, function(x) read_tsv(paste0(directory, x), col_names = "RSIDs"))

names(sumstat_rsids) = files

rsids = read_lines("derived_data/ANTE_RSIDS")

sumstat_rsids2 = lapply(sumstat_rsids, function(x) mutate(x, in_antegenes = RSIDs %in% rsids) %>% 
                          filter(!in_antegenes))

sumstats = read_tsv('derived_data/sumstats/PRACTICAL_PCa_Seibert2018.GRCh37.gz')
sumstats[sumstats$RSID %in% sumstat_rsids2$seibert_rsids$RSIDs, ]
