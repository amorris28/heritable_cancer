args = commandArgs(trailingOnly=TRUE)

model_file = args[1] # PHS file
inverse = args[2]
output_file = args[3] # output file name

# initialize packages
library(tidyverse)
library(lubridate)
library(readxl)

# constants
source("src/config.R")
source("src/functions.R")

# import data
avt_data = read_delim("/ess/p2521/data/durable/HuskMent/phenotypes/NPRdata/Tilleggsbestilling_Aug_2024/AVT.csv", delim = ";")
som_data = read_delim("/ess/p2521/data/durable/HuskMent/phenotypes/NPRdata/Tilleggsbestilling_Aug_2024/SOM.csv", delim = ";")
phv_data = read_delim("/ess/p2521/data/durable/HuskMent/phenotypes/NPRdata/Tilleggsbestilling_Aug_2024/PHV.csv", delim = ";")
death_data = read_excel("/ess/p2521/data/durable/HuskMent/phenotypes/DAAR_P210891_KOBLET.xlsx")
kreft_data = read_delim("/ess/p2521/data/durable/HuskMent/phenotypes/Kreftregisteret/Utlevert_Haavik/Utlevert_Haavik.csv", delim = ";")
age_dob_data = read_delim("/ess/p2521/data/durable/HuskMent/phenotypes/data_2023-09-28/HUSKMENT_data_withkeys_2023-09-28.csv", delim = ",")
birth_month_husk2 = read_delim("/ess/p2521/data/durable/HuskMent/phenotypes/HUSK_questionairre_data/WORK_QUERY_FOR_HUSK2_FMND_AAR.csv", delim = ",")

# Check last death dates
# no deaths after 12.2020
death_data %>% arrange(desc(DAAR), desc(DMND))

# IDs for each wave
husk1ids = age_dob_data$HUSKMentID
husk2ids = birth_month_husk2$HUSKMentID

#########################################################################
# Questionably motivated attempt #2
#########################################################################
## Step 1: All IDs
age_dob_data %>%
	select(huskmentid = HUSKMentID,
		   bfile_id = decode_subjectid) ->
		   id_key
sum(!unique(avt_data$lopenr) %in% id_key$huskmentid)
sum(!unique(som_data$lopenr) %in% id_key$huskmentid)
sum(!unique(phv_data$lopenr) %in% id_key$huskmentid)
sum(!unique(death_data$HUSKMENTID) %in% id_key$huskmentid)
sum(!unique(kreft_data$huskmentid) %in% id_key$huskmentid)
sum(!unique(birth_month_husk2$HUSKMentID) %in% id_key$huskmentid)
# It doesn't matter that these are missing because if I don't have DECODE IDs
# for them then I have no genotype data (or no way of connecting genotype data
# to them)

## Step 2: male IDs
# It looks like only cancer registry has gender data
kreft_data %>%
	select(huskmentid, kjoenn = KJOENN) %>%
	filter(kjoenn == "M") %>%
	pull(huskmentid) ->
	kreft_male_huskmentids

age_dob_data %>%
	select(huskmentid = HUSKMentID,
	       sex1 = Gender,
	       sex2 = kjona,
	       sex3 = kjonn) %>%
	mutate(sex1 = case_match(sex1, "M" ~ 1, "F" ~ 2))  %>%
	mutate(sex = coalesce(sex1, sex2, sex3)) %>%
	select(huskmentid, sex) %>%
    filter(sex == 1) %>%
    pull(huskmentid) ->
    male_huskmentids

all_male_huskmentids = unique(c(kreft_male_huskmentids, male_huskmentids))

#########################################################################
## Step 2: Age last obs
### Step 2a: birthdate

# ID DOB
age_dob_data %>%
	select(huskmentid = HUSKMentID, faar1, faar) %>%
	mutate(faar3 = if_else(is.na(faar1), faar, faar1)) %>%
    select(-faar, -faar1) %>%
	rename(faar1 = faar3) ->
	birth_year

birth_month_husk2 %>%
	rename(huskmentid = HUSKMentID, 
	       fmnd = invitasjonsregister_fmnd, 
	       faar2 = invitasjonsregister_faar) ->
	 birth_month

birth_year %>%
	full_join(birth_month, by = "huskmentid") ->
	dob_combined

# 10k subjects are missing from birth_month_husk2 and therefore have no birth months
# but they do have IDs from decode so they have been genotype
# This function will 
dob_combined %>%
	mutate(faar = coalesce(faar2, faar1)) %>%
	select(huskmentid, faar, fmnd) %>%
	mutate(faar = faar + 1900) %>%
	mutate(fmnd = if_else(is.na(fmnd), 1, fmnd)) %>%
	mutate(fdato = ym(paste(faar, fmnd))) %>%
	select(huskmentid, fdato) ->
	fodselsdato

## 2b: Date last obs
### registry data
avt_data %>%
	select(huskmentid = lopenr, dato = utdato) ->
	avt_data_proc
som_data %>%
	select(huskmentid = lopenr, dato = utdato) ->
	som_data_proc
phv_data %>%
	select(huskmentid = lopenr, dato = utdato) ->
	phv_data_proc

avt_data_proc %>%
	rbind(som_data_proc) %>%
	rbind(phv_data_proc) ->
	npr_dates

### death data
death_data %>%
	select(huskmentid = HUSKMENTID, daar = DAAR, dmnd = DMND) %>%
	mutate(dato = ym(paste(daar, dmnd, sep = "-"))) %>%
	select(huskmentid, dato) ->
	death_dates

### cancer data
kreft_data %>%
	select(huskmentid, diagnosedato = DIAGNOSEDATO) %>%
	mutate(dato = dmy(diagnosedato)) %>%
	select(huskmentid, dato) ->
	kreft_dates

npr_dates %>%
	rbind(death_dates) %>%
	rbind(kreft_dates) %>%
	group_by(huskmentid) %>%
	summarize(dato = max(dato)) %>%
	rename(date_lastobs = dato) ->
	censor_date

## 2c. Age last obs

fodselsdato %>%
	full_join(censor_date, by = "huskmentid") %>%
	mutate(age_lastobs = interval(fdato, date_lastobs) / years(1)) ->
	age_lastobs

#########################################################################
## Step 3: ID cases (and diagnosis date)
# Process avt_data
avt_data %>%
	select(huskmentid = lopenr, icd10, diagnose_dato = inndato) %>%
	mutate(case = grepl(icd10pca, icd10, fixed=TRUE)) %>%
	select(huskmentid, diagnose_dato, case) ->
	avt_data_cases

# Process som_data
som_data %>%
	select(huskmentid = lopenr, icd10, diagnose_dato = inndato) %>%
	mutate(case = grepl(icd10pca, icd10, fixed=TRUE)) %>%
	select(huskmentid, diagnose_dato, case) ->
	som_data_cases

# Process phv_data
phv_data %>%
	select(huskmentid = lopenr, icd10, diagnose_dato = inndato) %>%
	mutate(case = grepl(icd10pca, icd10, fixed=TRUE)) %>%
	select(huskmentid, diagnose_dato, case) ->
	phv_data_cases

# Process death data
death_data %>%
	select(huskmentid = HUSKMENTID, daar = DAAR, dmnd = DMND, alle_koder = ALLE_KODER) %>%
	mutate(diagnose_dato = ym(paste(daar, dmnd, sep = "-"))) %>%
	mutate(case = grepl(icd10pca, alle_koder, fixed=TRUE)) %>%
	select(huskmentid, diagnose_dato, case) ->
	death_data_cases

# Process kreft data
kreft_data %>%
	select(huskmentid, icd10 = ICD10, diagnose_dato = DIAGNOSEDATO) %>%
	mutate(case = grepl(icd10pca, icd10, fixed=TRUE)) %>%
	mutate(diagnose_dato = dmy(diagnose_dato)) %>%
	select(huskmentid, diagnose_dato, case) ->
	kreft_data_cases

avt_data_cases %>%
	rbind(som_data_cases) %>%
	rbind(phv_data_cases) %>%
	rbind(death_data_cases) %>%
	rbind(kreft_data_cases) %>%
	filter(case) %>%
	group_by(huskmentid) %>%
	mutate(first_diagnose_dato = min(diagnose_dato)) %>%
	ungroup() %>%
	filter(first_diagnose_dato == diagnose_dato) %>%
	select(-first_diagnose_dato) ->
	diagnose_date

age_lastobs %>%
	full_join(diagnose_date, by = "huskmentid") %>%
	mutate(status = if_else(is.na(case), 0, 1)) %>%
	mutate(age_diagnosis = interval(fdato, diagnose_dato) / years(1)) %>%
	mutate(age = if_else(status == 1, age_diagnosis, age_lastobs)) ->
	phenotype_data

#########################################################################
## Step 4: metastasis and age of metastasis
kreft_data %>%
	select(huskmentid, icd10 = ICD10, diagnose_dato = DIAGNOSEDATO, metastasis_code = METASTASE) %>%
	mutate(case = grepl(icd10pca, icd10, fixed=TRUE)) %>%
	mutate(metastasis_date = dmy(diagnose_dato)) %>%
	mutate(metastasis = metastasis_code %in% metastase_strict) %>%
	filter(case & metastasis) %>%
	select(huskmentid, metastasis_code, metastasis_date, metastasis) ->
	metastasis_data

phenotype_data %>%
	full_join(metastasis_data, by = "huskmentid") %>%
	mutate(age_of_metastasis = interval(fdato, diagnose_dato) / years(1)) %>%
	mutate(metastasis = if_else(is.na(metastasis), 0, metastasis)) %>%
	mutate(age_of_metastasis = if_else(metastasis == 1, age_of_metastasis, if_else(is.na(age_lastobs), age, age_lastobs))) %>%
	mutate(conversion_to_metastasis_years = age_of_metastasis - age_diagnosis) ->
	phenotype_data


# Filter for male subjects
phenotype_data %>%
	full_join(id_key, by = "huskmentid") %>%
	filter(huskmentid %in% all_male_huskmentids) ->
	phenotype_data

# Add PHS
phs_data = import_phs_data(model_file, inverse)

combined_data = process_combined_data(phenotype_data, phs_data)

# Save output
write_tsv(combined_data, output_file)
