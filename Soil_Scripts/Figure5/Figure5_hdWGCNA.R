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
library(hdWGCNA)
library(cowplot)
library(WGCNA)

current_dir <- trimws(getwd())
setwd(current_dir)
cat("Working directory:", current_dir, "\n")

parts <- strsplit(current_dir, .Platform$file.sep)[[1]]
project_root <- if (length(parts) >= 2 && all(tail(parts, 2) == c("Soil_Scripts", "Figure5"))) paste(parts[1:(length(parts)-2)], collapse = .Platform$file.sep) else current_dir
data_path <- file.path(project_root, "Soil_Matrix")
cat("Project root:", project_root, "\nData path:", data_path, "\n")

GLOBAL_PURITY_CUTOFF <- 0
name_list <- list.dirs(data_path, full.names = FALSE, recursive = FALSE)
seurat_list <- list()
initial_stats <- after_filter_stats <- after_seurat_stats <- data.frame(Sample = character(), Genes = integer(), Cells = integer())

pass_ttest <- function(x1, x2, x3) {
  values <- na.omit(as.numeric(c(x2, x3)))
  x1 <- as.numeric(x1)
  if (length(values) < 2 || var(values) == 0 || is.na(x1)) return(FALSE)
  x1 > max(values)
}

for (sample_name in name_list) {
  report_S_file <- file.path(data_path, paste0(sample_name, "_sc_taxonomy.report"))
  report_G_file <- file.path(data_path, paste0(sample_name, "_sc_taxonomy_G.report"))
  taxonomy_file <- report_S_file
  if (!all(file.exists(c(report_S_file, report_G_file, taxonomy_file)))) next
  
  bacteria_info <- read.delim(report_S_file, header = TRUE)
  colnames(bacteria_info) <- c("BC", "name", "taxonomy_id", "taxonomy_lvl", "reads", "all_reads")
  bacteria_info_G <- read.delim(report_G_file, header = TRUE)
  colnames(bacteria_info_G) <- c("BC", "name", "taxonomy_id", "taxonomy_lvl", "reads", "all_reads")
  tax_df <- read.delim(taxonomy_file, header = TRUE, stringsAsFactors = FALSE)
  
  expr_mat <- tryCatch(Read10X(file.path(data_path, sample_name)), error = function(e) NULL)
  if (is.null(expr_mat) || ncol(expr_mat) == 0) next
  
  initial_stats <- rbind(initial_stats, data.frame(Sample = sample_name, Genes = nrow(expr_mat), Cells = ncol(expr_mat)))
  
  filtered_tax_df <- tax_df[tax_df$fraction_total_reads >= GLOBAL_PURITY_CUTOFF, ]
  filtered_tax_df <- filtered_tax_df[apply(filtered_tax_df[, c("fraction_total_reads", "fraction_total_reads2", "fraction_total_reads3")], 1, \(x) pass_ttest(x[1], x[2], x[3])), ]
  keep_cells <- intersect(colnames(expr_mat), filtered_tax_df$barcode)
  expr_mat <- expr_mat[, keep_cells, drop = FALSE]
  
  after_filter_stats <- rbind(after_filter_stats, data.frame(Sample = sample_name, Genes = nrow(expr_mat), Cells = ncol(expr_mat)))
  if (ncol(expr_mat) == 0) next
  
  sample_obj <- CreateSeuratObject(expr_mat, min.cells = 10, min.features = 10)
  BC_info <- bacteria_info[match(colnames(sample_obj), bacteria_info$BC), ]
  BC_info_G <- bacteria_info_G[match(colnames(sample_obj), bacteria_info_G$BC), ]
  sample_obj$species_info <- BC_info$name
  sample_obj$genus_info <- BC_info_G$name
  sample_obj$sample <- sample_name
  sample_obj$crop <- factor(strsplit(sample_name, "_")[[1]][1])
  sample_obj$orig.ident <- sample_name
  sample_obj$gene_number <- Matrix::colSums(GetAssayData(sample_obj, slot = "counts") != 0)
  
  after_seurat_stats <- rbind(after_seurat_stats, data.frame(Sample = sample_name, Genes = nrow(sample_obj), Cells = ncol(sample_obj)))
  seurat_list[[sample_name]] <- sample_obj
}

write.csv(initial_stats, "initial_stats.csv", row.names = FALSE)
write.csv(after_filter_stats, "after_filter_stats.csv", row.names = FALSE)
write.csv(after_seurat_stats, "after_seurat_stats.csv", row.names = FALSE)

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

all_genes <- unique(unlist(lapply(seurat_list, rownames)))

for (nm in names(seurat_list)) {
  obj <- seurat_list[[nm]]
  missing <- setdiff(all_genes, rownames(obj))
  if (length(missing)) {
    counts_mat <- GetAssayData(obj, slot = "counts")
    zero <- Matrix(0, nrow = length(missing), ncol = ncol(obj), dimnames = list(missing, colnames(obj)))
    obj <- CreateSeuratObject(counts = rbind(counts_mat, zero), meta.data = obj@meta.data)
  }
  seurat_list[[nm]] <- obj
}

combined_seurat_object <- Reduce(merge, seurat_list)
cat("Merged object:", dim(combined_seurat_object), "\n")

write.csv(data.frame(Gene = rownames(combined_seurat_object)), "combined_seurat_all_genes.csv", row.names = FALSE)

combined_seurat_object <- JoinLayers(combined_seurat_object, assay = "RNA")
cat("RNA layers:", paste(Layers(combined_seurat_object[["RNA"]]), collapse = ", "), "\n")

combined_seurat_object <- subset(combined_seurat_object, subset = nCount_RNA > 0 & is.finite(nCount_RNA))

combined_seurat_object <- SCTransform(
  combined_seurat_object,
  vars.to.regress = NULL,
  verbose = FALSE,
  return.only.var.genes = FALSE,
  variable.features.n = 2000
)

combined_seurat_object <- RunPCA(combined_seurat_object, features = VariableFeatures(combined_seurat_object))
ElbowPlot(combined_seurat_object, ndims = 50, reduction = "pca")

find_pca_elbow <- function(seurat_obj, reduction = "pca", ndims = 50) {
  y <- seurat_obj[[reduction]]@stdev[1:ndims]
  x <- seq_along(y)
  x_norm <- (x - min(x)) / diff(range(x))
  y_norm <- (y - min(y)) / diff(range(y))
  elbow_point <- which.max(abs(y_norm - (1 - x_norm)) / sqrt(2))
  
  p <- ggplot(data.frame(PC = x, SD = y), aes(PC, SD)) +
    geom_point(size = 3, color = "#2E86AB") +
    geom_line(color = "#2E86AB") +
    geom_vline(xintercept = elbow_point, linetype = "dashed", color = "red") +
    annotate("text", x = elbow_point + 1, y = max(y), label = paste0("Elbow = PC", elbow_point), color = "red", size = 8, hjust = 0) +
    labs(title = "PCA Elbow Plot", x = "Principal Component", y = "Standard Deviation") +
    theme_bw(base_size = 20) +
    theme(
      plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 25, face = "bold"),
      axis.text = element_text(size = 25),
      panel.grid.major = element_line(color = "grey80", linetype = "dotted"),
      panel.grid.minor = element_blank()
    )
  
  print(p)
  cat("Suggested number of PCs:", elbow_point, "\n")
  elbow_point
}

elbow_pc <- find_pca_elbow(combined_seurat_object)

set.seed(1024)
combined_seurat_object <- combined_seurat_object %>%
  RunUMAP(reduction = "pca", dims = 1:elbow_pc) %>%
  FindNeighbors(reduction = "pca", dims = 1:elbow_pc) %>%
  FindClusters(resolution = 0.1)

saveRDS(combined_seurat_object, "combined_seurat_object.rds")
combined_seurat_object <- readRDS("combined_seurat_object.rds")

theme_set(theme_cowplot())
enableWGCNAThreads(nThreads = 16)

# Setup WGCNA
combined_seurat_object <- SetupForWGCNA(combined_seurat_object, gene_select = "fraction", fraction = 0.01, wgcna_name = "combined")
selected_genes <- GetWGCNAGenes(combined_seurat_object, wgcna_name = "combined")
cat("Number of selected genes:", length(selected_genes), "\n")

combined_seurat_object <- RunPCA(combined_seurat_object)
combined_seurat_object <- MetacellsByGroups(combined_seurat_object, reduction = "pca", k = 25, max_shared = 10, ident.group = "seurat_clusters")
combined_seurat_object <- NormalizeMetacells(combined_seurat_object)
combined_seurat_object <- subset(combined_seurat_object, features = selected_genes)

# Network construction
combined_seurat_object <- SetDatExpr(combined_seurat_object, assay = "SCT", layer = "data")
combined_seurat_object <- TestSoftPowers(combined_seurat_object, networkType = "signed")

plot_list <- PlotSoftPowers(combined_seurat_object)
pdf("FigureS_softpower.pdf", width = 10, height = 9)
print(wrap_plots(plot_list, ncol = 2))
dev.off()

power_table <- GetPowerTable(combined_seurat_object)
head(power_table)

combined_seurat_object <- ConstructNetwork(combined_seurat_object, tom_name = "combined_network", overwrite_tom = TRUE, minModuleSize = 30, mergeCutHeight = 0.25)

pdf("FigureS_dendrogram.pdf", width = 6, height = 4, useDingbats = FALSE)
PlotDendrogram(combined_seurat_object, main = "Gene dendrogram")
dev.off()

# Module analysis
combined_seurat_object <- ModuleEigengenes(combined_seurat_object)
combined_seurat_object <- ModuleConnectivity(combined_seurat_object)
PlotKMEs(combined_seurat_object, ncol = 3)

modules <- GetModules(combined_seurat_object) %>% filter(module != "grey")
module_df <- GetHubGenes(combined_seurat_object, n_hubs = 1000)
write.csv(module_df, "combined_module_genes.csv", row.names = FALSE, quote = FALSE)

# Module gene expression
module_gene_df <- data.frame(gene = unlist(modules$gene), module = rep(modules$module, sapply(modules$gene, length)))
module_genes <- intersect(module_gene_df$gene, rownames(GetAssayData(combined_seurat_object, assay = "SCT", slot = "data")))
expr_matrix <- GetAssayData(combined_seurat_object, assay = "SCT", slot = "data")
module_gene_expr_matrix <- expr_matrix[module_genes, , drop = FALSE]

crops <- combined_seurat_object$crop
crop_ids <- sort(unique(crops))
fraction_mat <- mean_mat <- matrix(0, nrow = length(module_genes), ncol = length(crop_ids), dimnames = list(module_genes, crop_ids))

for (tid in crop_ids) {
  cells <- names(crops[crops == tid])
  expr <- module_gene_expr_matrix[, cells, drop = FALSE]
  fraction_mat[, tid] <- rowSums(expr > 0) / length(cells)
  mean_mat[, tid] <- rowMeans(expr)
}

write.csv(fraction_mat, "module_gene_fraction_by_crop.csv", quote = FALSE)
write.csv(mean_mat, "module_gene_mean_expr_by_crop.csv", quote = FALSE)

# Module scores
combined_seurat_object <- ModuleExprScore(combined_seurat_object, n_genes = 50, method = "UCell")

plot_list <- ModuleFeaturePlot(combined_seurat_object, features = "scores", order = "shuffle", ucell = TRUE)
plot_list <- lapply(plot_list, \(p) p + geom_point(size = 0.01) + theme(plot.title = element_text(size = 30, face = "bold")))

pdf("Figure5a.pdf", width = 14, height = 10)
print(wrap_plots(plot_list, ncol = 3))
dev.off()

# Radar plots
combined_seurat_object$crop <- factor(combined_seurat_object$crop, levels = c("Soybean", "Rice", "Wheat"))

pdf("Figure5b.pdf", width = 12, height = 10)
ModuleRadarPlot(combined_seurat_object, group.by = "crop", axis.label.size = 5, grid.label.size = 0)
dev.off()

combined_seurat_object$seurat_clusters <- factor(combined_seurat_object$seurat_clusters)

pdf("Figure5c.pdf", width = 12, height = 10)
ModuleRadarPlot(combined_seurat_object, group.by = "seurat_clusters", axis.label.size = 5, grid.label.size = 0)
dev.off()

saveRDS(combined_seurat_object, "combined_seurat_object_hdWGCNA.rds")
combined_seurat_object <- readRDS("combined_seurat_object_hdWGCNA.rds")
cat("hdWGCNA analysis completed. Results saved to combined_seurat_object_hdWGCNA.rds\n")