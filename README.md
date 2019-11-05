# TCGA-LUSC
RNA-Seq analysis of TSGA-LUSC dataset using DESeq2

Script to download the source data from TCGA and create the raw counts matrix
*TCGA_LUSC_RNA_Seq_paired_analysis_01_source_data_download.R*

Script to clean the raw counts matrix (remove low rowSums) and assign gene names
*TCGA_LUSC_RNA_Seq_paired_analysis_02_raw_counts_matrix_cleaning.R*

Script to retrieve and manipulate demographic and tumour information from JSON file
*TCGA_LUSC_RNA_Seq_paired_analysis_03_demographics_from_JSON.R*

Script to create the dds object from counts matrix (countData), and demographics (colData) using DESeqDataSetFromMatrix
*TCGA_LUSC_RNA_Seq_paired_analysis_04_dds_object_from_matrix.R*

Script to analyse the vst-normalized data. Analysis include:
- PCA plots
- distance plot
*TCGA_LUSC_RNA_Seq_paired_analysis_05_vst_normalization_and_analysis.R*
