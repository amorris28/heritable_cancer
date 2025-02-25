source("src/00_config.R")

phsmodels = readr::read_lines(paste0(RAW_DATA, "phs-models.txt"))

##########################################
# Download and preprocess external sumstats files

files = tools::file_path_sans_ext(basename(phsmodels), compression = TRUE)

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

# New model from Roshan
PHS601_raw = paste0(EXTERNAL_DATA, "01_geneticScore.finalModel.snplist.xlsx")
PHS601_preprocfile = paste0(INTERIM_DATA, "PHS601.txt.gz")
PHS601_yamlfile = paste0(
                  tools::file_path_sans_ext(PHS601_preprocfile, compression = TRUE), 
                  ".yaml")

if (file.exists(PHS601_raw)) {
  PHS601 = readxl::read_excel(paste0(EXTERNAL_DATA, "01_geneticScore.finalModel.snplist.xlsx"))
  
  PHS601 = dplyr::select(PHS601,
                rsID = SNP,
                chr_name = Chromosome,
                chr_position = `Position (GRCh37)`,
                effect_allele = Effect,
                other_allele = Ref,
                effect_weight = Beta)
  
  
  if (!file.exists(PHS601_preprocfile)) {
    readr::write_tsv(PHS601, PHS601_preprocfile)
  }
  
  if (!file.exists(PHS601_yamlfile)) {
    readr::write_lines(paste0("cleansumstats_metafile_kind: minimal\ncol_SNP: rsID\ncol_CHR: chr_name\ncol_POS: chr_position\ncol_EffectAllele: effect_allele\ncol_OtherAllele: other_allele\ncol_BETA: effect_weight\npath_sumStats: ", basename(PHS601_preprocfile), "\nstats_Model: cox"),
                       PHS601_yamlfile)
  }
}
files = c(files, tools::file_path_sans_ext(basename(PHS601_preprocfile), compression = TRUE))
readr::write_lines(files, paste0(PROCESSED_DATA, "sumstats_basenames"))
