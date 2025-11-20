# Helper functions for single-cell analysis

#' Load and preprocess single-cell data
load_and_preprocess_data <- function() {
  # Sample information
  samples <- list(
    hiv_positive = c("ART1", "ART2", "ART3", "ART4", "ART5", 
                     "ART6", "ART7", "ART8", "ART9", "ART10"),
    healthy = c("HD1", "HD2", "HD3", "HD4", "HD5", "HD6", "HD7", "HD8", "HD9")
  )
  
  seurat_list <- list()
  
  # Load all samples
  for (sample_group in names(samples)) {
    for (sample in samples[[sample_group]]) {
      cat("Processing:", sample, "\n")
      
      # Load counts
      if (sample_group == "healthy") {
        data_path <- file.path(project_dir, "../0.public_healthy_tissue", 
                              paste0(sample, "_out/outs/filtered_feature_bc_matrix/"))
      } else {
        data_path <- file.path(project_dir, paste0(sample, "_out/outs/filtered_feature_bc_matrix/"))
      }
      
      counts <- Read10X(data_path)
      seurat_obj <- CreateSeuratObject(counts = counts, project = sample,
                                      min.cells = 3, min.features = 200)
      
      # Add mitochondrial percentage
      seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
      
      # Quality filtering
      qc_params <- get_qc_parameters(sample)
      seurat_obj <- subset(seurat_obj, 
                          subset = percent.mt < qc_params$mt_cutoff &
                          nFeature_RNA >= qc_params$min_features &
                          nFeature_RNA < qc_params$max_features)
      
      # Remove doublets if file exists
      seurat_obj <- remove_doublets(seurat_obj, sample)
      
      seurat_list[[sample]] <- seurat_obj
    }
  }
  
  return(seurat_list)
}

#' Get QC parameters for each sample
get_qc_parameters <- function(sample) {
  qc_params <- list(
    mt_cutoff = ifelse(grepl("HD", sample), 10, 50),
    min_features = ifelse(grepl("HD", sample), 1000, 200),
    max_features = ifelse(grepl("HD", sample), 7500, 9000)
  )
  return(qc_params)
}

#' Remove doublets if doublet file exists
remove_doublets <- function(seurat_obj, sample) {
  doublet_file <- file.path(project_dir, "1.results", paste0(sample, "_doublet.txt"))
  
  if (file.exists(doublet_file)) {
    doublet_data <- read.table(doublet_file, sep = ",", header = TRUE, row.names = 1)
    seurat_obj <- AddMetaData(seurat_obj, doublet_data)
    
    if ("predicted_doublets" %in% colnames(doublet_data)) {
      seurat_obj <- subset(seurat_obj, predicted_doublets == "False")
    } else if ("doublet_scores" %in% colnames(doublet_data)) {
      seurat_obj <- subset(seurat_obj, doublet_scores < 0.3)
    }
  }
  
  return(seurat_obj)
}

#' Integrate data using Harmony
integrate_data <- function(seurat_list) {
  # Merge all objects
  scRNA_merged <- merge(seurat_list[[1]], y = seurat_list[-1])
  
  # Standard preprocessing
  scRNA_merged <- NormalizeData(scRNA_merged, normalization.method = "LogNormalize", 
                               scale.factor = 10000)
  scRNA_merged <- FindVariableFeatures(scRNA_merged, selection.method = "vst", 
                                      nfeatures = 2000)
  scRNA_merged <- ScaleData(scRNA_merged, features = rownames(scRNA_merged))
  scRNA_merged <- RunPCA(scRNA_merged, features = VariableFeatures(scRNA_merged), 
                        verbose = FALSE)
  
  # Harmony integration
  scRNA_harmony <- RunHarmony(scRNA_merged, group.by.vars = "orig.ident")
  scRNA_harmony <- RunUMAP(scRNA_harmony, dims = 1:25, reduction = "harmony")
  
  return(scRNA_harmony)
}

#' Perform clustering analysis
cluster_analysis <- function(scRNA_harmony) {
  scRNA_harmony <- FindNeighbors(scRNA_harmony, dims = 1:25, reduction = "harmony")
  scRNA_harmony <- FindClusters(scRNA_harmony, resolution = 0.5)
  return(scRNA_harmony)
}

#' Annotate cell types
annotate_cell_types <- function(scRNA_harmony) {
  celltype_mapping <- list(
    '0' = 'CD4T cells',
    '3' = 'CD8T cells', 
    '5' = 'epithelial cells',
    '6' = 'epithelial cells',
    '7' = 'epithelial cells',
    '8' = 'epithelial cells',
    '12' = 'epithelial cells',
    '15' = 'endothelial cells',
    '9' = 'Stromal cells',
    '17' = 'Stromal cells', 
    '19' = 'Stromal cells',
    '1' = 'plasma cells',
    '4' = 'plasma cells',
    '2' = 'B cells',
    '13' = 'myeloid cells',
    '14' = 'myeloid cells',
    '10' = 'mast cells',
    '23' = 'mast cells',
    '11' = 'cycling cells',
    '16' = 'Tuft cells',
    '18' = 'ILCs',
    '20' = 'enteroendocrine L cells',
    '21' = 'iSMC',
    '22' = 'Neuron cells'
  )
  
  scRNA_harmony$celltype <- "Unknown"
  for (cluster_id in names(celltype_mapping)) {
    scRNA_harmony$celltype[scRNA_harmony$seurat_clusters == cluster_id] <- celltype_mapping[[cluster_id]]
  }
  
  # Add group information
  healthy_samples <- c("HD1", "HD2", "HD3", "HD4", "HD5", "HD6", "HD7", "HD8", "HD9")
  scRNA_harmony$group <- ifelse(scRNA_harmony$orig.ident %in% healthy_samples, "HD", "ART")
  
  return(scRNA_harmony)
}

#' Generate main figures
generate_main_figures <- function(scRNA_harmony) {
  # UMAP by cluster
  p1 <- DimPlot(scRNA_harmony, label = TRUE, repel = TRUE) +
    theme_minimal() +
    ggtitle("UMAP - Clusters")
  
  # UMAP by cell type
  colors <- c("#F19294", "#5D9BBE", "#F5B375", "#C0937E", "#67A59B",
              "#A4D38E", "#4A9D47", "#96C3D8", "#E45A5F", "#3477A9",
              "#BDA7CB", "#e49ac3", "#9983B7", "#CD9A99", "#DFCDE4",
              "#B61F8A", "#F58135", "#F7A96C", "#E78B75", "#DD9E82",
              "#DE4B3F", "#A65A34", "#E0C880", "#795FA3")
  
  p2 <- DimPlot(scRNA_harmony, group.by = "celltype", label = TRUE, cols = colors) +
    theme_minimal() +
    ggtitle("UMAP - Cell Types")
  
  # Save figures
  ggsave(file.path(figures_dir, "umap_clusters.pdf"), p1, width = 10, height = 8)
  ggsave(file.path(figures_dir, "umap_celltypes.pdf"), p2, width = 10, height = 8)
  
  # Cell proportion plot
  cell_prop_plot <- plot_cell_proportions(scRNA_harmony)
  ggsave(file.path(figures_dir, "cell_proportions.pdf"), cell_prop_plot, 
         width = 9, height = 5)
  
  # Marker gene heatmap
  marker_heatmap <- plot_marker_heatmap(scRNA_harmony)
  ggsave(file.path(figures_dir, "marker_heatmap.pdf"), marker_heatmap,
         width = 12, height = 10)
}

#' Plot cell proportions
plot_cell_proportions <- function(scRNA_harmony) {
  cell_ratio <- prop.table(table(scRNA_harmony$celltype, scRNA_harmony$orig.ident), margin = 2)
  cell_ratio_df <- as.data.frame(cell_ratio)
  
  # Order samples
  sample_order <- c("HD1", "HD2", "HD3", "HD4", "HD5", "HD6", 
                   "HD7", "HD8", "HD9", "ART1", "ART2", ART3", 
                   "ART4", "ART5", "ART6", "ART7", 
                   "ART8", "ART9", "ART10")
  
  cell_ratio_df$Var2 <- factor(cell_ratio_df$Var2, levels = sample_order)
  
  ggplot(cell_ratio_df) + 
    coord_flip() +
    geom_bar(aes(x = Var2, y = Freq, fill = Var1), 
             stat = "identity", width = 0.7, color = '#222222') + 
    theme_classic() +
    scale_fill_manual(values = colors) +
    labs(x = 'Sample', y = 'Proportion') +
    theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5))
}

#' Plot marker gene heatmap
plot_marker_heatmap <- function(scRNA_harmony) {
  marker_genes <- c("CD3E", "CD3D", "KRT8", "EPCAM", "JCHAIN", "MZB1", "XBP1", 
                   "BANK1", "MS4A1", "PDGFRB", "CXCL14", "CALD1", "PMP22", 
                   "CD14", "CD68", "MNDA", "KIT", "MS4A2", "MKI67", "TOP2A", 
                   "PLVAP", "PECAM1", "SH2D6", "IL7R", "CSF2", "PYY", "INSL5", 
                   "GCG", "ACTA2", "CDH19")
  
  DoHeatmap(scRNA_harmony, features = marker_genes, group.by = "celltype",
            assay = 'RNA', group.colors = colors) +
    scale_fill_gradientn(colors = c("white", "grey", "firebrick3"))
}
