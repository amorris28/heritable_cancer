source("src/config.R")

phsmodels = readr::read_lines(paste0(RAW_DATA, "phs-models.txt"))

##########################################
# Download and preprocess external sumstats files

for (url in phsmodels) {
  # import sumstats URLs and extract file names
  file = basename(url)
  destfile = paste0(EXTERNAL_DATA, file)
  yamlfile = paste0(INTERIM_DATA, 
                    tools::file_path_sans_ext(file, compression = TRUE), 
                    ".yaml")
  preprocfile = paste0(INTERIM_DATA, file)
  
  # Download files if they don't exist
  if (!file.exists(destfile)) {
    download.file(url, destfile)
    }
  
  # Remove leading comments
  if (!file.exists(preprocfile)) {
    temp = readr::read_tsv(destfile, comment = "#")
    temp = dplyr::select(temp, rsID, chr_name, chr_position, effect_allele, other_allele, effect_weight)
    readr::write_tsv(temp, preprocfile)
  }
  
  # Create YAML file
  if (!file.exists(yamlfile)) {
    readr::write_lines(paste0("cleansumstats_metafile_kind: minimal\ncol_SNP: rsID\ncol_CHR: chr_name\ncol_POS: chr_position\ncol_EffectAllele: effect_allele\ncol_OtherAllele: other_allele\ncol_BETA: effect_weight\npath_sumStats: ", file, "\nstats_Model: cox"),
                yamlfile)
    }
}
