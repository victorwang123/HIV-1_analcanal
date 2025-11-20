# Single Cell Sequencing Identifies HIV-specific Stem-like T Cells Associated with Immune Recovery During Suppressive Antiretroviral Therapy


This repository contains code for analyzing single-cell RNA sequencing data from HIV-infected and healthy tissue samples.

## Project Structure

| Path | Type | Description |
|------|------|-------------|
| `main_analysis.R` | File | Main analysis script |
| `config.R` | File | Configuration file |
| `README.md` | File | Project documentation |
| `R/` | Directory | R scripts and functions |
| `R/helper_functions.R` | File | Utility functions |
| `R/subset_analyses.R` | File | Subset analysis functions |
| `data/` | Directory | Input data (excluded from repo) |
| `results/` | Directory | Analysis results |
| `figures/` | Directory | Generated figures |


## Installation

1. Clone this repository:
```bash
git clone https://github.com/yourusername/hiv_tissue_scRNA.git
cd hiv_tissue_scRNA
install.packages(c("Seurat", "harmony", "dplyr", "ggplot2", "ggsci", "pheatmap"))
source("main_analysis.R")
