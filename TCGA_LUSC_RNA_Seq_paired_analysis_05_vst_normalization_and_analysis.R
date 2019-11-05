################################################################################
# Date: 2019-11-05
# Created by: RLADDACH
# Title: TCGA_LUSC_RNA_Seq_paired_analysis_05_vst_normalization_and_analysis.R
# Script to analyse the vst-normalized data. Analysis include:
# - PCA plots
# - distance plot
################################################################################

# loading libraries
library(DESeq2)       # Differential expression analysis
library(ggplot2)      # Visualisations
library(gridExtra)    # Visualisations -  multiple plots
library(pheatmap)     # Visualisations 
library(RColorBrewer) # Visualisations

# working directory setup
dir = "~/TCGA-LUSC-final"
setwd(dir)
rm(dir)

# loading the initial dds object
load("TCGA_LUSC_initial_dds_DESeqDataSetFromMatrix_paired_analysis.Rdata")

################################################################################
# vst transformation
vsd = vst(dds, blind=FALSE)

save(vsd, file = "TCGA_LUSC_vst_normalization_paired_analysis.Rdata")
write.csv(assay(vsd), "TCGA_LUSC_vst_normalization_paired_analysis.csv")

################################################################################
# PCA plots - one variable

# available intgroups:
# "condition", 
# "gender",
# "race",
# "smoking_status",
# "stage_com",
# "ajcc_path_T_com",
# "ajcc_path_N_com")

PCA_simple_condition = plotPCA(vsd, intgroup=c("condition")) + 
  ggtitle("PCA: condition") + 
  theme(legend.position="bottom", plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
PCA_simple_gender = plotPCA(vsd, intgroup=c("gender")) + 
  ggtitle("PCA: gender") + 
  theme(legend.position="bottom", plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
PCA_simple_race = plotPCA(vsd, intgroup=c("race")) + 
  ggtitle("PCA: race") + 
  theme(legend.position="bottom", plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
PCA_simple_smoking_status = plotPCA(vsd, intgroup=c("smoking_status")) + 
  ggtitle("PCA: smoking status") + 
  theme(legend.position="bottom", plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
PCA_simple_stage_com = plotPCA(vsd, intgroup=c("stage_com")) + 
  ggtitle("PCA: stage combined") + 
  theme(legend.position="bottom", plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
PCA_simple_ajcc_path_T_com = plotPCA(vsd, intgroup=c("ajcc_path_T_com")) + 
  ggtitle("PCA: AJCC path T combined") + 
  theme(legend.position="bottom", plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
PCA_simple_ajcc_path_N_com = plotPCA(vsd, intgroup=c("ajcc_path_N_com")) + 
  ggtitle("PCA: AJCC path N combined") + 
  theme(legend.position="bottom", plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

# combined plot of all 7 intgroups
grid.arrange(PCA_simple_condition,
             PCA_simple_gender,
             PCA_simple_race,
             PCA_simple_smoking_status,
             PCA_simple_stage_com,
             PCA_simple_ajcc_path_T_com,
             PCA_simple_ajcc_path_N_com,
             ncol=3)

# creating and saving plots of basic demographics
png("PCA_simple_combined_basic_demographics_paired_analysis.png", width = 1000, height = 750)
grid.arrange(PCA_simple_condition,
             PCA_simple_gender,
             PCA_simple_race,
             PCA_simple_smoking_status,
             ncol=2)
dev.off()

# creating and saving plots of staging information
png("PCA_simple_combined_staging_information_paired_analysis.png", width = 1000, height = 750)
grid.arrange(PCA_simple_condition,
             PCA_simple_stage_com,
             PCA_simple_ajcc_path_T_com,
             PCA_simple_ajcc_path_N_com,
             ncol=2)
dev.off()

################################################################################
# PCA plots - multi variable
# Am example of a complex graph with two intgroups - stage and race

PCA_complex_stage_vs_race = plotPCA(vsd, intgroup=c("stage_com", "race")) + 
  ggtitle("PCA: stage + race") +
  theme(legend.position="bottom", plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

################################################################################
# plots for the counts of stages
# combined stage
stage_com_plot = as.data.frame(t(table(vsd$stage_com)))
stage_com_plot = stage_com_plot[,2:3]
colnames(stage_com_plot) = c("combined_stage", "count")

png("counts_of_combined_stage_paired_analysis.png", width = 250, height = 200)
ggplot(data=stage_com_plot, aes(x=combined_stage, y=count)) + 
  geom_bar(stat="identity", width=0.7, fill="green4") +
  geom_text(aes(label=count), vjust=-0.1, color="black", size=4)
dev.off()

# AJCC pathologic T stage
ajcc_path_t_com_plot = as.data.frame(t(table(vsd$ajcc_path_T_com)))
ajcc_path_t_com_plot = ajcc_path_t_com_plot[,2:3]
colnames(ajcc_path_t_com_plot) = c("AJCC_T_stage", "count")

png("counts_of_AJCC_path_T_stage_paired_analysis.png", width = 250, height = 200)
ggplot(data=ajcc_path_t_com_plot, aes(x=AJCC_T_stage, y=count)) + 
  geom_bar(stat="identity", width=0.7, fill="green4") +
  geom_text(aes(label=count), vjust=-0.1, color="black", size=4)
dev.off()

# AJCC pathologic N stage
ajcc_path_n_com_plot = as.data.frame(t(table(vsd$ajcc_path_N_com)))
ajcc_path_n_com_plot = ajcc_path_n_com_plot[,2:3]
colnames(ajcc_path_n_com_plot) = c("AJCC_N_stage", "count")

png("counts_of_AJCC_path_N_stage_paired_analysis.png", width = 250, height = 200)
ggplot(data=ajcc_path_n_com_plot, aes(x=AJCC_N_stage, y=count)) + 
  geom_bar(stat="identity", width=0.7, fill="green4") +
  geom_text(aes(label=count), vjust=-0.1, color="black", size=3) +
  scale_x_discrete(limits=c("normal", "N0", "N1", "N2_3", "NX"))
dev.off()

################################################################################
# displaying sample distances
# calculating sample distance
sampleDists = dist(t(assay(vsd)))
?dist

# sampleDists
sampleDistMatrix = as.matrix(sampleDists)
rownames(sampleDistMatrix) = substr(paste(vsd$condition, colnames(vsd), sep="-"),1,15)
colnames(sampleDistMatrix) = NULL
colors = colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

png("eucledian_distance_from_vst_normalization_paired_analysis.png", width = 1300, height = 1300)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
dev.off()