###############################################################################
# Title: Untargeted Lipidomics: Data Analysis Workflow
# Author: Marya S.
# Date: 2023
#
# Description:
# This script performs data import, quality control, normalization,
# and downstream analyses of lipidomics datasets using the lipidR package.
# Analyses include:
#   - Principal Component Analysis (PCA)
#   - Hierarchical clustering
#   - Differential lipid abundance testing
#   - Lipid class enrichment analysis
#
# The pipeline is based on lipidR (v2.13.0; Mohamed et al., 2020) and
# LipidSig (Lin et al., 2021).
#
# Dependencies:
#   R (≥ 4.2.0)
#   lipidR (v2.13.0)
#   tidyverse
#   data.table
#   pheatmap
#   ggplot2
#
# References:
#   Mohamed, A.** *et al.* (2020). PMID: 32168452.
#   Lin, WJ.** *et al.* (2021). PMID: 34048582.
###############################################################################

# -----------------------------------------------------------------------------
# Installation and Setup
# -----------------------------------------------------------------------------

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("lipidr")

# If installation error occurs (macOS), link Fortran libraries in terminal:
# ln /usr/local/gfortran/lib/libgomp.1.dylib /Library/Frameworks/R.framework/Resources/lib/.
# ln /usr/local/gfortran/lib/libgomp.1.dylib /usr/local/lib/libgomp.1.dylib

install.packages("gmm")
library(lipidr)
use_interactive_graphics()

setwd("~/Desktop/NIHOxCam/Experiments/FSASD_Mice/Lipidomics/Analysis")

# -----------------------------------------------------------------------------
# Data Import and Quality Control (QC samples)
# -----------------------------------------------------------------------------

d1 <- as_lipidomics_experiment(
  read.csv("Untargeted Lipidomics Analysis Result_Positive Ion_CPLZ07222101-3_Edited1_QCOnly.csv")
)
d2 <- as_lipidomics_experiment(
  read.csv("Untargeted Lipidomics Analysis Result_Negative Ion_CPLZ07222101-3_Edited1_QCOnly.csv")
)
d1 <- add_sample_annotation(d1, "SampleAnnotation.csv")
d2 <- add_sample_annotation(d2, "SampleAnnotation.csv")

plot_samples(d1, type = "tic", log = TRUE)
plot_samples(d2, type = "tic", log = TRUE)
plot_molecules(d1, "cv", measure = "Area")
plot_molecules(d2, "cv", measure = "Area")
plot_lipidclass(d1, "boxplot")
plot_lipidclass(d2, "boxplot")

# -----------------------------------------------------------------------------
# Full Dataset (Experimental + QC samples)
# -----------------------------------------------------------------------------

d3 <- as_lipidomics_experiment(
  read.csv("Untargeted Lipidomics Analysis Result_Positive Ion_CPLZ07222101-3_Edited1.csv")
)
d4 <- as_lipidomics_experiment(
  read.csv("Untargeted Lipidomics Analysis Result_Negative Ion_CPLZ07222101-3_Edited1.csv")
)
d3 <- add_sample_annotation(d3, "SampleAnnotation.csv")
d4 <- add_sample_annotation(d4, "SampleAnnotation.csv")

plot_samples(d3, type = "tic", log = TRUE)
plot_samples(d4, type = "tic", log = TRUE)
plot_lipidclass(d3, "boxplot")
plot_lipidclass(d4, "boxplot")

# -----------------------------------------------------------------------------
# Data Summarization and Normalization
# -----------------------------------------------------------------------------

d_summarized1 <- summarize_transitions(d3, method = "average")
d_summarized2 <- summarize_transitions(d4, method = "average")

d_normalized1 <- normalize_pqn(d_summarized1, measure = "Area", exclude = "blank", log = TRUE)
d_normalized2 <- normalize_pqn(d_summarized2, measure = "Area", exclude = "blank", log = TRUE)

plot_samples(d_normalized1, "boxplot")
plot_samples(d_normalized2, "boxplot")

d_normalized_istd1 <- normalize_istd(d_summarized1, measure = "Area", exclude = "blank", log = TRUE)
d_normalized_istd2 <- normalize_istd(d_summarized2, measure = "Area", exclude = "blank", log = TRUE)

# -----------------------------------------------------------------------------
# Multivariate Analysis (PCA, OPLS-DA)
# -----------------------------------------------------------------------------

mvaresults1 <- mva(d_normalized1, measure = "Area", method = "PCA")
plot_mva(mvaresults1, color_by = "Group", components = c(1, 2))

mvaresults2 <- mva(d_normalized2, measure = "Area", method = "PCA")
plot_mva(mvaresults2, color_by = "Group", components = c(1, 2))

mvaresults1 <- mva(d_normalized1, method = "OPLS-DA", group_col = "Group", groups = c("WT", "R39CHOM"))
plot_mva(mvaresults1, color_by = "Group")

mvaresults2 <- mva(d_normalized2, method = "OPLS-DA", group_col = "Group", groups = c("WT", "R39CHOM"))
plot_mva(mvaresults2, color_by = "Group")

plot_mva_loadings(mvaresults1, color_by = "Class", top.n = 20)
plot_mva_loadings(mvaresults2, color_by = "Class", top.n = 20)

# -----------------------------------------------------------------------------
# Differential Analysis
# -----------------------------------------------------------------------------

library(limma)

de_results1 <- de_analysis(data = d_normalized1, R39CHOM - WT, measure = "Area")
de_results2 <- de_analysis(data = d_normalized2, R39CHOM - WT, measure = "Area")

head(de_results1)
head(de_results2)

significant_molecules(de_results1)
significant_molecules(de_results2)

plot_results_volcano(de_results1, show.labels = TRUE)
plot_results_volcano(de_results2, show.labels = TRUE)

# -----------------------------------------------------------------------------
# Enrichment Analysis
# -----------------------------------------------------------------------------

enrich_results1 <- lsea(de_results1, rank.by = "logFC")
significant_lipidsets(enrich_results1)
plot_enrichment(de_results1, significant_lipidsets(enrich_results1), annotation = "class")
plot_enrichment(de_results1, significant_lipidsets(enrich_results1), annotation = "unsat")

enrich_results2 <- lsea(de_results2, rank.by = "logFC")
significant_lipidsets(enrich_results2)
plot_enrichment(de_results2, significant_lipidsets(enrich_results2), annotation = "class")
plot_enrichment(de_results2, significant_lipidsets(enrich_results2), annotation = "unsat")

plot_chain_distribution(de_results1)
plot_chain_distribution(de_results2)

# -----------------------------------------------------------------------------
# Clustering Analysis
# -----------------------------------------------------------------------------

install.packages("ggdendro")
install.packages("iheatmapr_0.5.1.tar.gz", repos = NULL, type = "source")
library(iheatmapr)

plot_heatmap(d_normalized1, molecule_annotation = FALSE, sample_annotation = "Group2",
             cluster_rows = "hclust", cluster_cols = "hclust")

plot_heatmap(d_normalized2, molecule_annotation = FALSE, sample_annotation = "Group2",
             cluster_rows = "hclust", cluster_cols = "hclust")

###############################################################################
# End of Script
###############################################################################
