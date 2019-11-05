################################################################################
# Date: 2019-11-05
# Created by: RLADDACH
# Title: TCGA_LUSC_RNA_Seq_paired_analysis_04_dds_object_from_matrix.R
# Script to create the dds object from counts matrix (countData), 
# and demographics (colData) using DESeqDataSetFromMatrix
################################################################################

# loading libraries
library(DESeq2)

# working directory setup
dir = "~/TCGA-LUSC-final"
setwd(dir)
rm(dir)

# loading the demographics file
load("TCGA_LUSC_demographics_clean_from_JSON_paired_analysis.Rdata")

# loading the clean counts matrix
load("TCGA_LUSC_merged_full_counts_matrix_clean_paired_analysis.Rdata")

# creating a subset of information for the paired analysis using matched sample information
# matching demographics with samples (normal tissue then primary tumour)
demographics = data.frame(condition = c(rep('NT',49),rep('TP',49)))
rownames(demographics) = colnames(merged_full)
demographics = mutate(demographics, substr(colnames(merged_full), 1, 12))
colnames(demographics) = c("condition", "partial_id")

# merging the subset with the full list
demographics_combined = left_join(demographics,
                                  demographics_from_JSON,
                                  by=c("partial_id"="submitter_id")
)

# making the partial_id a factor for paired analysis
demographics_combined$partial_id = as.factor(demographics_combined$partial_id)

# subsetting values for tumour stages for normal tissue samples with "normal"
demographics_combined$ajcc_pathologic_n = as.character(demographics_combined$ajcc_pathologic_n)
demographics_combined$ajcc_pathologic_n[1:49] = "normal"
demographics_combined$ajcc_pathologic_n = as.factor(demographics_combined$ajcc_pathologic_n)

demographics_combined$ajcc_pathologic_m = as.character(demographics_combined$ajcc_pathologic_m)
demographics_combined$ajcc_pathologic_m[1:49] = "normal"
demographics_combined$ajcc_pathologic_m = as.factor(demographics_combined$ajcc_pathologic_m)

demographics_combined$ajcc_pathologic_t = as.character(demographics_combined$ajcc_pathologic_t)
demographics_combined$ajcc_pathologic_t[1:49] = "normal"
demographics_combined$ajcc_pathologic_t = as.factor(demographics_combined$ajcc_pathologic_t)

demographics_combined$stage_com = as.character(demographics_combined$stage_com)
demographics_combined$stage_com[1:49] = "normal"
demographics_combined$stage_com = as.factor(demographics_combined$stage_com)

demographics_combined$ajcc_path_t_com = as.character(demographics_combined$ajcc_path_t_com)
demographics_combined$ajcc_path_t_com[1:49] = "normal"
demographics_combined$ajcc_path_t_com = as.factor(demographics_combined$ajcc_path_t_com)

demographics_combined$ajcc_path_n_com = as.character(demographics_combined$ajcc_path_n_com)
demographics_combined$ajcc_path_n_com[1:49] = "normal"
demographics_combined$ajcc_path_n_com = as.factor(demographics_combined$ajcc_path_n_com)

save(demographics_combined, file="TCGA_LUSC_demographics_combined_98_subset_paired_analysis.Rdata")
write.csv(demographics_combined, "TCGA_LUSC_demographics_combined_98_subset_paired_analysis.csv")

###############################################################################
# working on the count data with DESeq2
# creating the countdata (counts matrix)
cts = as.matrix(merged_full)

# creating the coldata
# loading as character, however everything is converted to factor (required for paired analysis)
coldata = data.frame()
coldata = as.data.frame(cbind(as.character(demographics_combined$condition),
                              as.character(demographics_combined$partial_id),
                              as.character(demographics_combined$gender),
                              as.character(demographics_combined$race),
                              as.character(demographics_combined$smoking_status),
                              as.character(demographics_combined$stage_com),
                              as.character(demographics_combined$ajcc_path_t_com),
                              as.character(demographics_combined$ajcc_path_n_com))
                        )

colnames(coldata) = c("condition", 
                      "partial_id",
                      "gender",
                      "race",
                      "smoking_status",
                      "stage_com",
                      "ajcc_path_T_com",
                      "ajcc_path_N_com")

rownames(coldata) = colnames(cts)

# sanity check
str(coldata)

# creating the dds object
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ partial_id + condition)
dds

save(dds, file="TCGA_LUSC_initial_dds_DESeqDataSetFromMatrix_paired_analysis.Rdata")