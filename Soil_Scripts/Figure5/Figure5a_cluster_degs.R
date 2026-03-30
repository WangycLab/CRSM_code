# ===============================================
# Author : ZHENG XINGHAI
# Date   : 2025-09-30
# Project: Soil mscRNA-seq Analysis
# ===============================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(harmony)
library(tidyverse)
library(tidyr)
library(ggrepel)
library(glmGamPoi)
library(RColorBrewer)
library(ggsci)
library(viridis)
library(Matrix)
library(monocle3)
library(SeuratWrappers)
library(patchwork)
library(irlba)

# Set working directory to the current location
current_dir <- trimws(getwd())
setwd(current_dir)
cat("Working directory set to:", current_dir, "\n")

# Determine project root
parts <- strsplit(current_dir, .Platform$file.sep)[[1]]
if (length(parts) >= 2 && all(tail(parts, 2) == c("Soil_Scripts", "Figure5"))) {
  project_root <- paste(parts[1:(length(parts)-2)], collapse = .Platform$file.sep)
} else {
  project_root <- current_dir
}
cat("Project root:", project_root, "\n")

# Define data directory relative to project root
data_path <- file.path(project_root, "Soil_Matrix")
cat("Data path set to:", data_path, "\n")

# Set global purity cutoff for filtering cells
GLOBAL_PURITY_CUTOFF <- 0

# Initialize containers for samples, statistics, and expression matrices
name_list <- list.dirs(data_path, full.names = FALSE, recursive = FALSE)
seurat_list <- list()
initial_stats <- data.frame(Sample = character(), Genes = integer(), Cells = integer())
after_filter_stats <- data.frame(Sample = character(), Genes = integer(), Cells = integer())
after_seurat_stats <- data.frame(Sample = character(), Genes = integer(), Cells = integer())
expr_count_raw <- list()
expr_count_filtered <- list()
gene_number_list <- list()

# Define t-test based function to filter low-purity cells
pass_ttest <- function(fraction_total_reads, fraction_total_reads2, fraction_total_reads3) {
  values <- as.numeric(c(fraction_total_reads2, fraction_total_reads3))
  fraction_total_reads <- as.numeric(fraction_total_reads)
  values <- values[!is.na(values)]
  if (length(values) < 2 || var(values) == 0 || is.na(fraction_total_reads)) {
    return(FALSE)
  }
  test <- t.test(values, mu = fraction_total_reads)
  test$p.value <= 1 && fraction_total_reads > max(values)
}

# Process each sample independently and apply filtering
for (sample_name in name_list) {
  report_S_file <- file.path(data_path, paste0(sample_name, "_sc_taxonomy.report"))
  report_G_file <- file.path(data_path, paste0(sample_name, "_sc_taxonomy_G.report"))
  taxonomy_file <- file.path(data_path, paste0(sample_name, "_sc_taxonomy.report"))
  
  if (!file.exists(report_S_file) || !file.exists(report_G_file) || !file.exists(taxonomy_file)) next
  
  # Load species and genus information
  bacteria_info <- read.delim(report_S_file, header = TRUE)
  colnames(bacteria_info) <- c("BC","name","taxonomy_id","taxonomy_lvl","reads","all_reads")
  bacteria_info_G <- read.delim(report_G_file, header = TRUE)
  colnames(bacteria_info_G) <- c("BC","name","taxonomy_id","taxonomy_lvl","reads","all_reads")
  tax_df <- read.delim(taxonomy_file, header = TRUE, stringsAsFactors = FALSE)
  
  # Load expression matrix
  sample_dir <- file.path(data_path, sample_name)
  expr_mat <- tryCatch(Read10X(sample_dir), error = function(e) NULL)
  if (is.null(expr_mat) || ncol(expr_mat) == 0) next
  
  # Record initial expression matrix statistics
  initial_stats <- rbind(
    initial_stats, 
    data.frame(Sample = sample_name, Genes = nrow(expr_mat), Cells = ncol(expr_mat))
  )
  
  # Filter cells based on purity and t-test
  filtered_tax_df <- tax_df[tax_df$fraction_total_reads >= GLOBAL_PURITY_CUTOFF, ]
  filtered_tax_df <- filtered_tax_df[apply(filtered_tax_df[, c("fraction_total_reads", "fraction_total_reads2", "fraction_total_reads3")], 1, function(x) pass_ttest(x[1], x[2], x[3])), ]
  keep_cells <- intersect(colnames(expr_mat), filtered_tax_df$barcode)
  expr_mat <- expr_mat[, keep_cells, drop = FALSE]
  
  # Record post-filtering statistics
  after_filter_stats <- rbind(
    after_filter_stats, 
    data.frame(Sample = sample_name, Genes = nrow(expr_mat), Cells = ncol(expr_mat))
  )
  if (ncol(expr_mat) == 0) next
  
  # Create Seurat object and annotate metadata
  sample_obj <- CreateSeuratObject(expr_mat, min.cells = 10, min.features = 10)
  BC_info <- bacteria_info[match(colnames(sample_obj), bacteria_info$BC), ]
  BC_info_G <- bacteria_info_G[match(colnames(sample_obj), bacteria_info_G$BC), ]
  sample_obj$species_info <- BC_info$name
  sample_obj$genus_info <- BC_info_G$name
  sample_obj$sample <- sample_name
  parts <- strsplit(sample_name, "_")[[1]]
  sample_obj$crop <- factor(parts[1])
  sample_obj$orig.ident <- sample_name
  sample_obj$gene_number <- Matrix::colSums(GetAssayData(sample_obj, slot = "counts") != 0)
  
  # Record Seurat object statistics
  after_seurat_stats <- rbind(
    after_seurat_stats, 
    data.frame(Sample = sample_name, Genes = nrow(sample_obj), Cells = ncol(sample_obj))
  )
  
  seurat_list[[sample_name]] <- sample_obj
}

# Export per-sample statistics for reference
write.csv(initial_stats, "initial_stats.csv", row.names = FALSE)
write.csv(after_filter_stats, "after_filter_stats.csv", row.names = FALSE)
write.csv(after_seurat_stats, "after_seurat_stats.csv", row.names = FALSE)

# Add sample name as prefix to cell barcodes if missing
for (sample_name in names(seurat_list)) {
  obj <- seurat_list[[sample_name]]
  bc <- colnames(obj)
  need_prefix <- !grepl(paste0("^", sample_name, "_"), bc)
  if (any(need_prefix)) {
    bc[need_prefix] <- paste0(sample_name, "_", bc[need_prefix])
    colnames(obj) <- bc
    rownames(obj@meta.data) <- bc
  }
  seurat_list[[sample_name]] <- obj
}

# Harmonize gene sets across samples and merge all Seurat objects
all_genes <- unique(unlist(lapply(seurat_list, rownames)))
for (nm in names(seurat_list)) {
  obj <- seurat_list[[nm]]
  missing <- setdiff(all_genes, rownames(obj))
  if (length(missing) > 0) {
    zero <- matrix(0, nrow = length(missing), ncol = ncol(obj), dimnames = list(missing, colnames(obj)))
    counts_mat <- GetAssayData(obj, slot = "counts")
    obj <- CreateSeuratObject(counts = rbind(counts_mat, zero), meta.data = obj@meta.data)
  }
  seurat_list[[nm]] <- obj
}
combined_seurat_object <- Reduce(merge, seurat_list)
cat("All samples merged, object dimension:", dim(combined_seurat_object), "\n")

# Export all gene names from merged object
all_genes_combined <- rownames(combined_seurat_object)
genes_df <- data.frame(Gene = all_genes_combined)
write.csv(genes_df, file = "combined_seurat_all_genes.csv", row.names = FALSE)
cat("Total genes exported:", length(all_genes_combined), "\n")

# Remove cells with zero or invalid RNA counts
combined_seurat_object <- subset(combined_seurat_object, subset = nCount_RNA > 0 & is.finite(nCount_RNA))

# Normalize expression and perform PCA
combined_seurat_object <- SCTransform(combined_seurat_object, vars.to.regress = NULL, verbose = FALSE, return.only.var.genes = FALSE, variable.features.n = 2000)
combined_seurat_object <- RunPCA(combined_seurat_object, features = VariableFeatures(combined_seurat_object))

# Visualize PCA using an elbow plot
ElbowPlot(combined_seurat_object, ndims = 50, reduction = "pca")

# Function to automatically detect the PCA elbow point
find_pca_elbow <- function(seurat_obj, reduction = "pca", ndims = 50, axis_title_size = 25, axis_text_size = 25, axis_font = "sans", title_size = 30, title_font = "sans") {
  if (!reduction %in% names(seurat_obj@reductions)) stop("Specified reduction not found")
  pca_sd <- seurat_obj[[reduction]]@stdev[1:ndims]
  
  x <- 1:ndims
  y <- pca_sd
  x_norm <- (x - min(x)) / (max(x) - min(x))
  y_norm <- (y - min(y)) / (max(y) - min(y))
  
  distance <- abs(y_norm - (1 - x_norm)) / sqrt(2)
  elbow_point <- which.max(distance)
  
  plot_df <- data.frame(PC = x, SD = y)
  p <- ggplot(plot_df, aes(x = PC, y = SD)) +
    geom_point(size = 3, color = "#2E86AB") +
    geom_line(color = "#2E86AB") +
    geom_vline(xintercept = elbow_point, linetype = "dashed", color = "red") +
    annotate("text", x = elbow_point + 1, y = max(y), label = paste0("Elbow = PC", elbow_point), color = "red", size = 8, hjust = 0) +
    labs(title = "PCA Elbow Plot", x = "Principal Component", y = "Standard Deviation") +
    theme_bw(base_size = 20) +
    theme(
      plot.title = element_text(family = title_font, size = title_size, face = "bold", hjust = 0.5),
      axis.title = element_text(family = axis_font, size = axis_title_size, face = "bold"),
      axis.text = element_text(family = axis_font, size = axis_text_size),
      panel.grid.major = element_line(color = "grey80", linetype = "dotted"),
      panel.grid.minor = element_blank()
    )
  print(p)
  cat("Suggested number of PCs:", elbow_point, "\n")
  return(elbow_point)
}

# Detect PCA elbow point automatically
elbow_pc <- find_pca_elbow(combined_seurat_object)

# Perform clustering and UMAP embedding based on selected PCs
set.seed(1024)
combined_seurat_object <- combined_seurat_object %>%
  RunUMAP(reduction = "pca", dims = 1:elbow_pc) %>%
  FindNeighbors(reduction = "pca", dims = 1:elbow_pc) %>%
  FindClusters(resolution = 0.1)

# Export final cell metadata
Idents(combined_seurat_object) <- "seurat_clusters"
md <- combined_seurat_object@meta.data
output_df <- data.frame(
  barcode = rownames(md),
  cluster = as.character(md$seurat_clusters),
  sample = md$orig.ident,
  crop = md$crop,
  species = md$species_info,
  genus = md$genus_info,
  gene_number = md$gene_number,
  stringsAsFactors = FALSE
)
write.table(output_df, file = "cell_metadata.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Perform differential expression analysis across clusters
combined_seurat_object <- PrepSCTFindMarkers(combined_seurat_object)

markers <- FindAllMarkers(
  combined_seurat_object,
  assay = "SCT",
  only.pos = TRUE,
  min.pct = 0.05,
  logfc.threshold = 0.25
)

markers <- markers %>%
  arrange(cluster, desc(avg_log2FC))

write.table(
  markers,
  file = "all_clusters_DEGs.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Extract normalized expression matrix and cluster identities
counts <- GetAssayData(combined_seurat_object, assay = "SCT", slot = "data")
clusters <- combined_seurat_object$seurat_clusters

markers$gene <- gsub("_", "-", markers$gene)
genes <- unique(markers$gene)
genes <- genes[genes %in% rownames(counts)]

cluster_ids <- sort(unique(clusters))
cluster_names <- paste0("Cluster_", cluster_ids)

fraction_mat <- matrix(0, nrow = length(genes), ncol = length(cluster_ids),
                       dimnames = list(genes, cluster_names))
mean_mat <- matrix(0, nrow = length(genes), ncol = length(cluster_ids),
                   dimnames = list(genes, cluster_names))

# Calculate expression fraction and mean expression per cluster
for (cid in cluster_ids) {
  cluster_cells <- names(clusters[clusters == cid])
  sub_counts <- counts[genes, cluster_cells, drop = FALSE]
  
  fraction <- rowSums(sub_counts > 0) / length(cluster_cells)
  mean_expr <- rowMeans(sub_counts)
  
  fraction_mat[, paste0("Cluster_", cid)] <- fraction
  mean_mat[, paste0("Cluster_", cid)] <- mean_expr
}

write.csv(fraction_mat, "gene_fraction_by_cluster.csv", quote = FALSE)
write.csv(mean_mat, "gene_mean_expr_by_cluster.csv", quote = FALSE)