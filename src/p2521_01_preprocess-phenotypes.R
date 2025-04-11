# initialize packages
library(tidyverse)
library(lubridate)

# constants
icd10pca = "C61"
source("src/00_config.R")

# import data
avt_data = read_delim("/ess/p2521/data/durable/HuskMent/phenotypes/NPRdata/Tilleggsbestilling_Aug_2024/AVT.csv", delim = ";")
som_data = read_delim("/ess/p2521/data/durable/HuskMent/phenotypes/NPRdata/Tilleggsbestilling_Aug_2024/SOM.csv", delim = ";")
phv_data = read_delim("/ess/p2521/data/durable/HuskMent/phenotypes/NPRdata/Tilleggsbestilling_Aug_2024/PHV.csv", delim = ";")
death_data = read_delim("/ess/p2521/cluster/users/p2521-andrewmo/data/DAAR_P210891_KOBLET.csv", delim = ",")
kreft_data = read_delim("/ess/p2521/data/durable/HuskMent/phenotypes/Kreftregisteret/Utlevert_Haavik/Utlevert_Haavik.csv", delim = ";")
age_dob_data = read_delim("/ess/p2521/data/durable/HuskMent/phenotypes/data_2023-09-28/HUSKMENT_data_withkeys_2023-09-28.csv", delim = ",")
birth_month_husk2 = read_delim("/ess/p2521/data/durable/HuskMent/phenotypes/HUSK_questionairre_data/WORK_QUERY_FOR_HUSK2_FMND_AAR.csv", delim = ",")

# Check last death dates
# no deaths after 12.2020
death_data %>% arrange(desc(DAAR), desc(DMND))

# functions
my_union = function(x, y) {
	if (is.na(x)) {
		if (is.na(y)) {
			return(NA)
		} else {
			return(y)
			}
	} else {
		 return(x)
	}
}

# ID Key
age_dob_data %>%
	select(subjectid = HUSKMentID,
	       decode_subjectid) ->
	id_key


###############################################################################################
all_ids = c(avt_data$lopenr, som_data$lopenr, phv_data$lopenr, death_data$HUSKMENTID, kreft_data$huskmentid) %>%
	unique()

husk1ids = age_dob_data$HUSKMentID
husk2ids = birth_month_husk2$HUSKMentID

##################################################################################################
# ID sex

age_dob_data %>%
	select(subjectid = HUSKMentID,
	       sex1 = Gender,
	       sex2 = kjona,
	       sex3 = kjonn) %>%
	mutate(sex1 = case_match(sex1, "M" ~ 1, "F" ~ 2)) ->
	gender_data_proc

gender_data_proc$sex <- sapply(1:nrow(gender_data_proc), function(x) my_union(gender_data_proc[[x,"sex1"]], gender_data_proc[[x,"sex2"]]))
gender_data_proc$sex <- sapply(1:nrow(gender_data_proc), function(x) my_union(gender_data_proc[[x,"sex"]], gender_data_proc[[x,"sex3"]]))

gender_data_proc %>%
	select(subjectid, sex) ->
	gender_data_proc

gender_data_proc %>%
        filter(sex == 1) %>%
        pull(subjectid) ->
        male_ids

#######################################################################################
# Pre-process phenotype data

# Process avt_data
avt_data %>%
	select(subjectid = lopenr, icd10, inndato, utdato) %>%
	mutate(case = icd10 == icd10pca) ->
	avt_data_proc

# Process som_data
som_data %>%
	select(subjectid = lopenr, icd10, inndato, utdato) %>%
	mutate(case = icd10 == icd10pca) ->
	som_data_proc

# Process phv_data
phv_data %>%
	select(subjectid = lopenr, icd10, inndato, utdato) %>%
	mutate(case = icd10 == icd10pca) ->
	phv_data_proc

# Process death data
death_data %>%
	select(subjectid = HUSKMENTID, daar = DAAR, dmnd = DMND, alle_koder = ALLE_KODER) %>%
	mutate(case = grepl(icd10pca, alle_koder, fixed=TRUE)) ->
	death_data_proc

# Process kreft data
kreft_data %>%
	select(subjectid = huskmentid, icd10 = ICD10, diagnosedato = DIAGNOSEDATO) %>%
	mutate(case = icd10 == icd10pca) ->
	kreft_data_proc

###################################################################################################
# ID DOB
age_dob_data %>%
	select(subjectid = HUSKMentID, faar1, faar) %>%
	mutate(faar3 = if_else(is.na(faar1), faar, faar1)) %>%
      	select(-faar, -faar1) %>%
	rename(faar1 = faar3) ->
	age_dob_data_proc

birth_month_husk2 %>%
	rename(subjectid = HUSKMentID, 
	       fmnd = invitasjonsregister_fmnd, 
	       faar2 = invitasjonsregister_faar) ->
	 birth_month_husk2_proc

age_dob_data_proc %>%
	full_join(birth_month_husk2_proc, by = "subjectid") ->
	dob_combined

dob_combined$faar <- sapply(1:nrow(dob_combined), function(x) my_union(dob_combined[[x,"faar1"]], dob_combined[[x,"faar2"]]))

dob_combined %>%
	select(subjectid, faar, fmnd) %>%
	mutate(faar = faar + 1900) %>%
	mutate(fmnd = if_else(is.na(fmnd), 1, fmnd)) %>%
	mutate(fdato = ym(paste(faar, fmnd))) %>%
	select(subjectid, fdato) ->
	fodselsdato

# Identify all cases
avt_data_proc %>%
	filter(case) %>%
	pull(subjectid) %>%
	unique() ->
	avt_data_pca_cases

som_data_proc %>%
	filter(case) %>%
	pull(subjectid) %>%
	unique() ->
	som_data_pca_cases

phv_data_proc %>%
	filter(case) %>%
	pull(subjectid) %>%
	unique() ->
	phv_data_pca_cases

death_data_proc %>%
	filter(case) %>%
	pull(subjectid) ->
	death_data_pca_cases

kreft_data_proc %>%
	filter(case) %>%
	pull(subjectid) %>%
	unique() ->
	kreft_data_pca_cases

c(avt_data_pca_cases, 
  som_data_pca_cases, 
  phv_data_pca_cases, 
  death_data_pca_cases,
  kreft_data_pca_cases) %>%
	unique() ->
	combined_cases

#####################################################################################
# Cases
# Keep all visits with ICD10 code for PCA.
# Get the earliest date (inndato, diagnosedato, death date)
avt_data_proc %>%
	rbind(som_data_proc) %>%
	rbind(phv_data_proc) %>%
	filter(case) %>%
	select(subjectid, diagnosedato = inndato) ->
	npr_data_proc_cases

death_data_proc %>%
	filter(case) %>%
	mutate(diagnosedato = ym(paste(daar, dmnd, sep = "-"))) %>%
	select(subjectid, diagnosedato)	->
	death_data_proc_cases

kreft_data_proc %>%
	filter(case) %>%
	mutate(diagnosedato = dmy(diagnosedato)) %>%
	select(subjectid, diagnosedato) ->
	kreft_data_proc_cases

npr_data_proc_cases %>%
	rbind(death_data_proc_cases) %>%
	rbind(kreft_data_proc_cases) ->
	all_data_proc_cases

all_data_proc_cases %>%
	group_by(subjectid) %>%
	summarize(dato = min(diagnosedato)) %>%
	mutate(status = 1) ->
	diagnosedato


# Controls
# Here, I remove subjects whose subject ID is in the list of cases
# Then, get the latest date for each control subject, this will be their censoring date
# Either last date observed (utdato/diagnosedato) or death date

# Modify this to pull the last date observed for all cases AND controls then filter by controls
avt_data_proc %>%
	rbind(som_data_proc) %>%
	rbind(phv_data_proc) %>%
	# filter(!subjectid %in% combined_cases) %>%
	select(subjectid, dato = utdato) ->
	npr_dates
	# npr_data_controls

death_data_proc %>%
	# filter(!subjectid %in% combined_cases) %>%
	mutate(dato = ym(paste(daar, dmnd, sep = "-"))) %>%
	select(subjectid, dato) ->
	death_dates
	# death_data_proc_controls

kreft_data_proc %>%
	# filter(!subjectid %in% combined_cases) %>%
	mutate(dato = dmy(diagnosedato)) %>%
	select(subjectid, dato) ->
	kreft_dates
	# kreft_data_proc_controls

# npr_data_controls %>%
# 	rbind(death_data_proc_controls) %>%
# 	rbind(kreft_data_proc_controls) ->
# 	all_data_proc_controls

npr_dates %>%
	rbind(death_dates) %>%
	rbind(kreft_dates) ->
	all_dates
	# all_data_proc_controls

# all_data_proc_controls %>%
all_dates %>%
	group_by(subjectid) %>%
	summarize(dato = max(dato)) %>%
	# mutate(status = 0) %>%
	filter(subjectid %in% male_ids) ->
	sensurdato
write_tsv(sensurdato, paste0(PROCESSED_DATA, "husk_sensurdato.tsv"))

sensurdato %>%
	filter(!subjectid %in% combined_cases) %>%
	mutate(status = 0) ->
	control_sensurdato

#################################################################################
diagnosedato %>%
	rbind(control_sensurdato) %>%
	full_join(fodselsdato, by = join_by(subjectid)) %>%
	mutate(age = interval(fdato, dato) / years(1)) %>%
	select(subjectid, age, status) %>%
	na.omit() ->
	phenotype_data

phenotype_data %>%
	left_join(id_key, by = join_by(subjectid)) %>%
	mutate(subjectid = decode_subjectid) %>%
	select(-decode_subjectid) ->
	phenotype_data

write_tsv(phenotype_data, paste0(PROCESSED_DATA, "husk_phenotypes.tsv"))


tibble(X1 = 0, X2 = phenotype_data$subjectid) ->
	pfile_ids

# Write list of subject ids to processed folder
# 
write_tsv(pfile_ids, paste0(INTERIM_DATA, "p2521_bfile_ids.txt"), col_names = FALSE)

# Export death data
death_data_proc %>%
	left_join(id_key, by = "subjectid") %>%
	mutate(subjectid = decode_subjectid) %>%
	select(-decode_subjectid) ->
	death_data_proc

write_tsv(death_data_proc, paste0(PROCESSED_DATA, "death_data_proc.tsv"))

write_tsv(diagnosedato, paste0(PROCESSED_DATA, "diagnosedato.tsv"))

write_tsv(id_key, paste0(PROCESSED_DATA, "id_key.tsv"))



