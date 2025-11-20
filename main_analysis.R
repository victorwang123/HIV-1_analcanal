#!/usr/bin/env Rscript
# HIV Tissue Single-Cell Analysis
# Clean version for GitHub

# Load required packages
required_packages <- c("Seurat", "harmony", "SeuratDisk", "ggsci", "dplyr", 
                      "ggplot2", "pheatmap", "viridis", "ggpubr", "Nebulosa")
invisible(lapply(required_packages, function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}))

# Set directories
project_dir <- "tissue_out/"
results_dir <- file.path(project_dir, "2.new_results")
figures_dir <- file.path(project_dir, "figures")

dir.create(results_dir, showWarnings = FALSE)
dir.create(figures_dir, showWarnings = FALSE)

setwd(project_dir)

# Source helper functions
source("R/helper_functions.R")

# Main analysis function
main <- function() {
  cat("Starting HIV tissue single-cell analysis...\n")
  
  # Load and preprocess data
  cat("1. Loading and preprocessing data...\n")
  seurat_list <- load_and_preprocess_data()
  
  # Merge and integrate data
  cat("2. Data integration with Harmony...\n")
  scRNA_harmony <- integrate_data(seurat_list)
  
  # Cluster analysis
  cat("3. Clustering analysis...\n")
  scRNA_harmony <- cluster_analysis(scRNA_harmony)
  
  # Cell type annotation
  cat("4. Cell type annotation...\n")
  scRNA_harmony <- annotate_cell_types(scRNA_harmony)
  
  # Save results
  cat("5. Saving results...\n")
  saveRDS(scRNA_harmony, file.path(results_dir, "integrated_seurat_object.rds"))
  
  # Generate main figures
  cat("6. Generating figures...\n")
  generate_main_figures(scRNA_harmony)
  
  cat("Analysis complete!\n")
  return(scRNA_harmony)
}

# Run main analysis
if (!file.exists(file.path(results_dir, "integrated_seurat_object.rds"))) {
  scRNA_harmony <- main()
} else {
  cat("Loading pre-computed results...\n")
  scRNA_harmony <- readRDS(file.path(results_dir, "integrated_seurat_object.rds"))
}

# Subset analyses
if (FALSE) {  # Set to TRUE to run subset analyses
  source("R/subset_analyses.R")
  
  # CD4 T cell analysis
  cd4_subset <- subset_and_analyze(scRNA_harmony, clusters = 0, 
                                   name = "CD4_Tcells", resolution = 1.5)
  
  # CD8 T cell analysis  
  cd8_subset <- subset_and_analyze(scRNA_harmony, clusters = 3,
                                   name = "CD8_Tcells", resolution = 2)
  
  # Epithelial cell analysis
  epithelial_subset <- subset_and_analyze(scRNA_harmony, clusters = c(5,6,7,8,12),
                                         name = "epithelial", resolution = 0.5)
}
