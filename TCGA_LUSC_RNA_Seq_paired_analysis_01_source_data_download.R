################################################################################
# Date: 2019-11-04
# Created by: RLADDACH
# Title: TCGA_LUSC_RNA_Seq_paired_analysis_01_source_data_download.R
# Script to retrieve data for TCGA-LUSC analysis and create a raw counts matrix 
################################################################################

# loading libraries
library(TCGAbiolinks) # Downloading the source data from TCGA

# working directory setup
dir = "~/TCGA-LUSC-final"
setwd(dir)
rm(dir)

# obtaining the initial information
query = GDCquery(project = "TCGA-LUSC", 
                 data.category = "Transcriptome Profiling", 
                 data.type = "Gene Expression Quantification", 
                 workflow.type = "HTSeq - Counts")

samplesDown = getResults(query,cols = c("cases"))

# matching cases and control for paired analysis (NT = normal,TP = tumour)
matched <- TCGAquery_MatchedCoupledSampleTypes(samplesDown,c("NT","TP"))

# subsetting NT and TP
dataSmNT <- TCGAquery_SampleTypes(barcode = matched,
                                  typesample = "NT")
dataSmTP <- TCGAquery_SampleTypes(barcode = matched,
                                  typesample = "TP")

# downloading the subsets
queryDown <- GDCquery(project = "TCGA-LUSC", 
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "HTSeq - Counts", 
                      barcode = c(dataSmNT, dataSmTP))
GDCdownload(queryDown, method = "api", files.per.chunk = 10)

# reading and preparing the RangedSummarizedExperimemt object
dataPrep_ranged <- GDCprepare(query = queryDown, 
                              save = TRUE, 
                              save.filename = "TCGA_LUSC_HTSeq_Counts_full_paired_analysis.Rdata")

# Raw matrix count preparation + AAIC plot
dataPrep_matrix <- TCGAanalyze_Preprocessing(object = dataPrep_ranged, 
                                             cor.cut = 0.6,
                                             datatype = "HTSeq - Counts",
                                             filename = "TCGA_LUSC_AAI_correlation_low_res_paired_analysis.png",
                                             width = 1500,
                                             height = 1500)

write.csv(dataPrep_matrix, "TCGA_LUSC_dataPrep_raw_matrix_full_paired_analysis.csv")
save(dataPrep_matrix, file = "TCGA_LUSC_dataPrep_raw_matrix_full_paired_analysis.Rdata")
