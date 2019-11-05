################################################################################
# Date: 2019-11-05
# Created by: RLADDACH
# Title: TCGA_LUSC_RNA_Seq_paired_analysis_03_demographics_from_JSON.R
# Script to retrieve and manipulate demographic and tumour information from JSON file
# Source JSON file from:
# https://portal.gdc.cancer.gov/exploration?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.project_id%22%2C%22value%22%3A%5B%22TCGA-LUSC%22%5D%7D%7D%5D%7D
# Select Clinical --> JSON
################################################################################

# loading libraries
library(jsonlite)
library(dplyr)

# working directory setup
dir = "~/TCGA-LUSC-final"
setwd(dir)

# checking the name of the json file
list.files(path = dir, pattern = "json")
rm(dir)

# reading JSON file to R
JSON_demo = fromJSON("clinical.cases_selection.2019-11-05.json")

################################################################################
# retrieving data from diagnoses
# year_of_diagnosis
get_year_of_diagnosis = function(x){
  if ('year_of_diagnosis' %in% names(x)){
    return(as.character(x['year_of_diagnosis']))
  }
  else{
    return(NA)
  }
}
year_of_diagnosis_from_JSON = as.character(lapply(JSON_demo$diagnoses,get_year_of_diagnosis))

# primary_diagnosis
get_primary_diagnosis = function(x){
  if ('primary_diagnosis' %in% names(x)){
    return(as.character(x['primary_diagnosis']))
  }
  else{
    return(NA)
  }
}
primary_diagnosis_from_JSON = as.character(lapply(JSON_demo$diagnoses,get_primary_diagnosis))

# tumor_stage
get_tumor_stage = function(x){
  if ('tumor_stage' %in% names(x)){
    return(as.character(x['tumor_stage']))
  }
  else{
    return(NA)
  }
}
tumor_stage_from_JSON = as.character(lapply(JSON_demo$diagnoses,get_tumor_stage))

# age_at_diagnosis
get_age_at_diagnosis = function(x){
  if ('age_at_diagnosis' %in% names(x)){
    return(as.character(x['age_at_diagnosis']))
  }
  else{
    return(NA)
  }
}
age_at_diagnosis_from_JSON = as.character(lapply(JSON_demo$diagnoses,get_age_at_diagnosis))

# ajcc_pathologic_n
get_ajcc_pathologic_n = function(x){
  if ('ajcc_pathologic_n' %in% names(x)){
    return(as.character(x['ajcc_pathologic_n']))
  }
  else{
    return(NA)
  }
}
ajcc_pathologic_n_from_JSON = as.character(lapply(JSON_demo$diagnoses,get_ajcc_pathologic_n))

# ajcc_pathologic_m
get_ajcc_pathologic_m = function(x){
  if ('ajcc_pathologic_m' %in% names(x)){
    return(as.character(x['ajcc_pathologic_m']))
  }
  else{
    return(NA)
  }
}
ajcc_pathologic_m_from_JSON = as.character(lapply(JSON_demo$diagnoses,get_ajcc_pathologic_m))

# ajcc_pathologic_t
get_ajcc_pathologic_t = function(x){
  if ('ajcc_pathologic_t' %in% names(x)){
    return(as.character(x['ajcc_pathologic_t']))
  }
  else{
    return(NA)
  }
}
ajcc_pathologic_t_from_JSON = as.character(lapply(JSON_demo$diagnoses,get_ajcc_pathologic_t))

# ajcc_staging_system_edition
get_ajcc_staging_system_edition = function(x){
  if ('ajcc_staging_system_edition' %in% names(x)){
    return(as.character(x['ajcc_staging_system_edition']))
  }
  else{
    return(NA)
  }
}
ajcc_staging_system_edition_from_JSON = as.character(lapply(JSON_demo$diagnoses,get_ajcc_staging_system_edition))

# ajcc_pathologic_stage
get_ajcc_pathologic_stage = function(x){
  if ('ajcc_pathologic_stage' %in% names(x)){
    return(as.character(x['ajcc_pathologic_stage']))
  }
  else{
    return(NA)
  }
}
ajcc_pathologic_stage_from_JSON = as.character(lapply(JSON_demo$diagnoses,get_ajcc_pathologic_stage))

################################################################################
# retrieving data from exposures
# years_smoked
get_years_smoked = function(x){
  if ('years_smoked' %in% names(x)){
    return(as.character(x['years_smoked']))
  }
  else{
    return(NA)
  }
}
years_smoked_from_JSON = as.character(lapply(JSON_demo$exposures,get_years_smoked))

# pack_years_smoked
get_pack_years_smoked = function(x){
  if ('pack_years_smoked' %in% names(x)){
    return(as.character(x['pack_years_smoked']))
  }
  else{
    return(NA)
  }
}
pack_years_smoked_from_JSON = as.character(lapply(JSON_demo$exposures,get_pack_years_smoked))

################################################################################
# data from demographics can be directly retrieved

################################################################################
# creating a combined data frame

demographics_from_JSON = data.frame(submitter_id = JSON_demo$demographic$submitter_id,
                                    gender = JSON_demo$demographic$gender,
                                    race = JSON_demo$demographic$race,
                                    ethnicity = JSON_demo$demographic$ethnicity,
                                    year_of_birth = as.numeric(JSON_demo$demographic$year_of_birth),
                                    age_at_index = as.numeric(JSON_demo$demographic$age_at_index),
                                    vital_status = JSON_demo$demographic$vital_status,
                                    year_of_diagnosis = as.numeric(year_of_diagnosis_from_JSON),
                                    primary_diagnosis = primary_diagnosis_from_JSON,
                                    tumor_stage = tumor_stage_from_JSON,
                                    age_at_diagnosis = as.numeric(age_at_diagnosis_from_JSON),
                                    ajcc_pathologic_n = ajcc_pathologic_n_from_JSON,
                                    ajcc_pathologic_m = ajcc_pathologic_m_from_JSON,
                                    ajcc_pathologic_t = ajcc_pathologic_t_from_JSON,
                                    ajcc_staging_system_edition = ajcc_staging_system_edition_from_JSON,
                                    years_smoked = as.numeric(years_smoked_from_JSON),
                                    pack_years_smoked = as.numeric(pack_years_smoked_from_JSON)
)

################################################################################
# adding derived columns

# smoking_status from pack_years_smoked
demographics_from_JSON = demographics_from_JSON %>%
  mutate(smoking_status = case_when(
    is.na(pack_years_smoked) ~ "unknown",
    pack_years_smoked > 0 ~ "smoker"
  ))

# collapsing stage_com into 3 main stages
# stage III and IV combined due to low numbers in stage IV (n=7)
demographics_from_JSON = demographics_from_JSON %>%
  mutate(stage_com = case_when(
    tumor_stage == "stage ia" ~ "stage_I",
    tumor_stage == "stage ib" ~ "stage_I",
    tumor_stage == "stage ii" ~ "stage_II",
    tumor_stage == "stage iia" ~ "stage_II",
    tumor_stage == "stage iib" ~ "stage_II",
    tumor_stage == "stage iiia" ~ "stage_III_IV",
    tumor_stage == "stage iiib" ~ "stage_III_IV",
    tumor_stage == "stage iv" ~ "stage_III_IV"
  ))

# collapsing ajcc_path_t_com into 4 main stages
demographics_from_JSON = demographics_from_JSON %>%
  mutate(ajcc_path_t_com = case_when(
    ajcc_pathologic_t == "stage ia" ~ "T1",
    ajcc_pathologic_t == "T1" ~ "T1",
    ajcc_pathologic_t == "T1a" ~ "T1",
    ajcc_pathologic_t == "T1b" ~ "T1",
    ajcc_pathologic_t == "T2" ~ "T2",
    ajcc_pathologic_t == "T2a" ~ "T2",
    ajcc_pathologic_t == "T2b" ~ "T2",
    ajcc_pathologic_t == "T3" ~ "T3",
    ajcc_pathologic_t == "T4" ~ "T4"
  ))

# collapsing ajcc_path_n_com into 5 groups due to low numbers in N3 (n=5)
demographics_from_JSON = demographics_from_JSON %>%
  mutate(ajcc_path_n_com = case_when(
    ajcc_pathologic_n == "N0" ~ "N0",
    ajcc_pathologic_n == "N1" ~ "N1",
    ajcc_pathologic_n == "N2" ~ "N2_3",
    ajcc_pathologic_n == "N3" ~ "N2_3",
    ajcc_pathologic_n == "NX" ~ "NX"
  ))

################################################################################
# Cleaning columns - removing extra text ("_demographics") or empty spaces

demographics_from_JSON$submitter_id = gsub("_demographic", "", demographics_from_JSON$submitter_id)
demographics_from_JSON$race = gsub(" ", "_", demographics_from_JSON$race)
demographics_from_JSON$ethnicity = gsub(" ", "_", demographics_from_JSON$ethnicity)
demographics_from_JSON$primary_diagnosis = gsub(",", "_", demographics_from_JSON$primary_diagnosis)

save(demographics_from_JSON, file="TCGA_LUSC_demographics_clean_from_JSON_paired_analysis.Rdata")
write.csv(demographics_from_JSON, "TCGA_LUSC_demographics_clean_from_JSON_paired_analysis.csv")