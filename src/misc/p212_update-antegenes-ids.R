library(tidyverse)

source("src/00_config.R")

fam = read_tsv('/cluster/projects/p212/antenor_pipeline/antegenes_imputed_allchr-temporary.fam',
	col_names = c("FID", "IID", "FA", "MO", "SEX", "PHE"))

fam %>%
	mutate(IID = sub("(.*?_.*?)_.*", "\\1", IID)) ->
	new_fam

bind_cols(fam[1:2], new_fam[1:2]) ->
	output

write_tsv(output, paste0(INTERIM_DATA, "update-ids-antegenes.tsv"), col_names = FALSE)

############################
# Output
# 1. Old family ID
# 2. Old within-family ID
# 3. New family ID
# 4. New within-family ID


