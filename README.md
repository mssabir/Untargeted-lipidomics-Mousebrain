# Untargeted Lipidomics: Data Analysis Workflow

This repository contains an R-based pipeline (using LipidR) for the processing and statistical analysis of untargeted lipidomics data from mouse brain tissues. 
The workflow was designed for studies investigating lipid dysregulation in Free Sialic Acid Storage Disorder (FSASD) mouse models.

---

## Overview

This workflow performs:

1. **Data import and quality control (QC)**
2. **Normalization**
   - Probabilistic Quotient Normalization (PQN)
   - Internal standard normalization (optional)
3. **Exploratory data analysis**
   - Principal Component Analysis (PCA)
   - Hierarchical clustering
4. **Differential lipid abundance testing**
   - Limma-based univariate comparison
5. **Lipid class enrichment and pathway analysis**
   - Lipid class and unsaturation-based enrichment
   - Pathway mapping via LipidSig (web interface only)

---

## Dependencies

| Package | Version | Reference |
|----------|----------|------------|
| R | ≥ 4.2.0 | [R Project](https://www.r-project.org/) |
| lipidR | 2.13.0 | Mohamed *et al.*, 2020 |
| tidyverse | ≥ 2.0 | — |
| data.table | ≥ 1.14 | — |
| pheatmap | ≥ 1.0 | — |
| ggplot2 | ≥ 3.5 | — |
| limma | ≥ 3.60 | — |
| ggdendro | ≥ 0.1 | — |
| iheatmapr | ≥ 0.5.1 | — |

---

## Installation

### 1. Install dependencies

```r
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("lipidr")

install.packages(c("tidyverse", "data.table", "pheatmap", "ggplot2", "limma", "ggdendro"))
install.packages("iheatmapr_0.5.1.tar.gz", repos = NULL, type = "source")
```

If an error occurs during `lipidr` installation on macOS, link the Fortran library in Terminal:

```bash
sudo ln /usr/local/gfortran/lib/libgomp.1.dylib /Library/Frameworks/R.framework/Resources/lib/.
sudo ln /usr/local/gfortran/lib/libgomp.1.dylib /usr/local/lib/libgomp.1.dylib
```

---

## Directory Structure

```
├── Lipidomics_Analysis_Workflow.R      # Main analysis script
├── SampleAnnotation.csv                # Sample metadata file
├── SampleAnnotation2.csv               # Alternative annotation (three-group comparison)
├── Untargeted Lipidomics Analysis Result_Positive Ion_*.csv
├── Untargeted Lipidomics Analysis Result_Negative Ion_*.csv
├── README.md                           # Documentation
```

---

## Running the Analysis

1. Open **RStudio** and set your working directory:
   ```r
   setwd("~/Desktop/NIHOxCam/Experiments/FSASD_Mice/Lipidomics/Analysis")
   ```
2. Run the workflow:
   ```r
   source("Lipidomics_Analysis_Workflow.R")
   ```

3. The script will:
   - Import the raw `.csv` lipidomics data files  
   - Add sample annotations  
   - Perform PQN normalization  
   - Conduct PCA and OPLS-DA  
   - Execute differential and enrichment analyses  
   - Generate volcano plots, enrichment maps, and heatmaps  

---

## Input File Requirements

| File | Description |
|------|--------------|
| `Untargeted Lipidomics Analysis Result_*.csv` | Raw output from mass spectrometry (positive and negative ion modes). |
| `SampleAnnotation.csv` | Metadata file defining sample names, groupings, and conditions. Must include `SampleID` and `Group` columns. |
| `SampleAnnotation2.csv` | Optional extended annotation for three-group or region-specific analyses. |

---

## Output

The workflow generates the following output files and figures:

- PCA and OPLS-DA score plots  
- Volcano plots (differential abundance)
- Lipid class enrichment bar plots  
- Hierarchical clustering heatmaps  
- Normalized lipid intensity matrices  

Figures are automatically displayed in RStudio and can be saved manually via the RStudio “Plots” panel.

---

## Citation

If you use this pipeline, please cite:

- **Mohamed, A.** *et al.* (2020). *lipidR: A software package for lipidomics data analysis.* *Journal of Proteome Research*, **19**(7), 2890–2898.  
- **Lin, K.** *et al.* (2021). *LipidSig: Comprehensive lipidomic data analysis and visualization.* *Nucleic Acids Research*, **49**(W1), W336–W344.  

---
