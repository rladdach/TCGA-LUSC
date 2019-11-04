################################################################################
# Date: 2019-11-04
# Created by: RLADDACH
# Title: TCGA_LUSC_RNA_Seq_paired_analysis_02_raw_counts_matrix_cleaning.R
# Script to clean the raw counts matrix (remove low rowSums) and assign gene names
################################################################################

# loading libraries
library(AnnotationDbi)
library(org.Hs.eg.db)

# working directory setup
dir = "~/TCGA-LUSC-final"
setwd(dir)
rm(dir)

# loading the raw counts matrix 
load("TCGA_LUSC_dataPrep_raw_matrix_full_paired_analysis.Rdata")
dataPrep_matrix = as.data.frame(dataPrep_matrix)

# initial size of the matrix = 56512 x 98
dim(dataPrep_matrix) 

# removing rows with rowSums<10 = 11629 rows removed
dataPrep_matrix = dataPrep_matrix[rowSums(dataPrep_matrix)>=10, ]

# data frame with row names from dataPrep_matrix
ensembl_gene = data.frame(ensembl_gene_id=rownames(dataPrep_matrix),
                          stringsAsFactors=FALSE)

# assigning gene names (SYMBOL) from AnnotationDbi/org.Hs.eg.db
ensembl_gene$symbol <- mapIds(org.Hs.eg.db,
                     keys=ensembl_gene$ensembl_gene_id,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

# adding a column for merging
dataPrep_matrix$ensembl_gene_id = rownames(dataPrep_matrix)

# merging gene symbols with raw count matrix
merged_full = merge(ensembl_gene, dataPrep_matrix, by="ensembl_gene_id")

# removing NA values
merged_full = merged_full[complete.cases(merged_full),]
dim(merged_full) # 23206 x 100 --> 21677 rows removed

# saving the mapping as multiple ensemble codes are mapped to one gene
ensemble_to_gene_mapping = merged_full[,c(1,2)]
save(ensemble_to_gene_mapping, file="TCGA_LUSC_dataPrep_ensemble_to_gene_mapping_paired_analysis.Rdata")

# cleaning - assigning row names and removing unwanted columns
merged_full$symbol =  make.unique(merged_full$symbol) # 47 genes duplicated
rownames(merged_full) = merged_full$symbol
merged_full$ensembl_gene_id = NULL
merged_full$symbol = NULL

# final size of the raw counts matrix 23206 x 98
dim(merged_full) 

save(merged_full, file = "TCGA_LUSC_merged_full_counts_matrix_clean_paired_analysis.Rdata")
write.csv(merged_full, "TCGA_LUSC_merged_full_counts_matrix_clean_paired_analysis.csv")