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
library(edgeR)

current_dir <- trimws(getwd())
setwd(current_dir)
cat("Working directory:", current_dir, "\n")

parts <- strsplit(current_dir, .Platform$file.sep)[[1]]
project_root <- if (length(parts) >= 2 && all(tail(parts, 2) == c("Soil_Scripts", "Figure6"))) paste(parts[1:(length(parts)-2)], collapse = .Platform$file.sep) else current_dir
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

set.seed(1024)

target_species <- "Burkholderia cepacia"

if (!"species_info" %in% colnames(combined_seurat_object@meta.data))
  stop("species_info column not found in meta.data. Please check upstream annotation steps.")

species_cells <- rownames(combined_seurat_object@meta.data)[combined_seurat_object@meta.data$species_info == target_species]
if (!length(species_cells)) stop(paste("No cells found for species:", target_species))

combined_seurat_object <- subset(combined_seurat_object, cells = species_cells)
cat("Remaining cells after filtering", target_species, ":", ncol(combined_seurat_object), "\n")

combined_seurat_object <- RunPCA(combined_seurat_object, features = VariableFeatures(combined_seurat_object))

find_pca_elbow <- function(seurat_obj, reduction = "pca", ndims = 50, axis_title_size = 25, axis_text_size = 25,
                           axis_font = "sans", title_size = 30, title_font = "sans") {
  if (!reduction %in% names(seurat_obj@reductions)) stop("Specified reduction not found in the Seurat object.")
  pca_sd <- seurat_obj[[reduction]]@stdev[seq_len(ndims)]
  x <- seq_len(ndims); y <- pca_sd
  x_norm <- (x - min(x)) / diff(range(x)); y_norm <- (y - min(y)) / diff(range(y))
  elbow_point <- which.max(abs(y_norm - (1 - x_norm)) / sqrt(2))
  p <- ggplot(data.frame(PC = x, SD = y), aes(PC, SD)) +
    geom_point(size = 3, color = "#2E86AB") +
    geom_line(color = "#2E86AB") +
    geom_vline(xintercept = elbow_point, linetype = "dashed", color = "red") +
    annotate("text", x = elbow_point + 1, y = max(y), label = paste0("Elbow = PC", elbow_point),
             color = "red", size = 8, hjust = 0) +
    labs(title = "PCA Elbow Plot", x = "Principal Component", y = "Standard Deviation") +
    theme_bw(base_size = 20) +
    theme(
      plot.title = element_text(family = title_font, size = title_size, face = "bold", hjust = 0.5),
      axis.title = element_text(family = axis_font, size = axis_title_size, face = "bold"),
      axis.text = element_text(family = axis_font, size = axis_text_size),
      panel.grid.major = element_line(color = "grey80", linetype = "dotted"),
      panel.grid.minor = element_blank()
    )
  print(p); cat("Suggested number of principal components:", elbow_point, "\n")
  elbow_point
}

elbow_pc <- find_pca_elbow(combined_seurat_object)

combined_seurat_object <- combined_seurat_object %>%
  RunUMAP(reduction = "pca", dims = 1:elbow_pc) %>%
  FindNeighbors(reduction = "pca", dims = 1:elbow_pc) %>%
  FindClusters(resolution = 0.1)

cat("UMAP and clustering completed. Final object dimension:", dim(combined_seurat_object), "\n")

DefaultAssay(combined_seurat_object) <- "SCT"
Idents(combined_seurat_object) <- "crop"
combined_seurat_object <- PrepSCTFindMarkers(combined_seurat_object)

markers_by_crop <- FindAllMarkers(combined_seurat_object, assay = "SCT", only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.25)
write.table(markers_by_crop, "DEGs_by_crop_species.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

counts <- GetAssayData(combined_seurat_object, assay = "SCT", slot = "data")
crops <- setNames(combined_seurat_object$crop, Cells(combined_seurat_object))

markers_by_crop$gene <- gsub("_", "-", markers_by_crop$gene)
genes <- intersect(unique(markers_by_crop$gene), rownames(counts))
cat("Number of DEG genes used for crop statistics:", length(genes), "\n")

crop_ids <- sort(unique(crops))
fraction_mat_crop <- mean_mat_crop <- matrix(0, nrow = length(genes), ncol = length(crop_ids),
                                             dimnames = list(genes, crop_ids))

for (st in crop_ids) {
  crop_cells <- names(crops)[crops == st]
  sub_counts <- counts[genes, crop_cells, drop = FALSE]
  fraction_mat_crop[, st] <- rowSums(sub_counts > 0) / length(crop_cells)
  mean_mat_crop[, st] <- rowMeans(sub_counts)
}

write.csv(fraction_mat_crop, "gene_fraction_by_crop.csv", quote = FALSE)
write.csv(mean_mat_crop, "gene_mean_expr_by_crop.csv", quote = FALSE)

Idents(combined_seurat_object) <- "seurat_clusters"
md_species <- combined_seurat_object@meta.data
if (!"seurat_clusters" %in% colnames(md_species))
  stop("seurat_clusters not found. Please ensure FindClusters() has been run.")

output_df_species <- data.frame(
  barcode = rownames(md_species),
  cluster = as.character(md_species$seurat_clusters),
  sample = md_species$orig.ident,
  crop = md_species$crop,
  stringsAsFactors = FALSE
)

write.table(output_df_species, "cell_metadata_species.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

clusters <- levels(Idents(combined_seurat_object))
cluster_colors <- setNames(colorRampPalette(brewer.pal(12, "Paired"))(length(clusters)), clusters)

umap_coords <- as.data.frame(Embeddings(combined_seurat_object, "umap"))
x_range <- diff(range(umap_coords$umap_1)); y_range <- diff(range(umap_coords$umap_2))
x_min <- min(umap_coords$umap_1); x_max <- max(umap_coords$umap_1)
y_min <- min(umap_coords$umap_2); y_max <- max(umap_coords$umap_2)

axis_len_x <- 0.2 * x_range; axis_len_y <- 0.2 * y_range
axis_x0 <- x_min - 0.05 * x_range; axis_y0 <- y_min - 0.05 * y_range
axis_arrow <- arrow(type = "closed", length = unit(5, "mm"))

p1 <- DimPlot(combined_seurat_object, reduction = "umap", label = FALSE, cols = cluster_colors, pt.size = 1, raster = FALSE) +
  labs(title = "UMAP by Cluster") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 40, family = "Arial"),
    axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
    legend.text = element_text(size = 20, family = "Arial"),
    legend.key.height = unit(2, "lines"), panel.grid = element_blank(),
    plot.margin = margin(10, 10, 90, 90, "pt")
  ) +
  coord_cartesian(clip = "off") +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 10))) +
  annotate("segment", x = axis_x0, xend = axis_x0 + axis_len_x, y = axis_y0, yend = axis_y0, linewidth = 2, arrow = axis_arrow) +
  annotate("segment", x = axis_x0, xend = axis_x0, y = axis_y0, yend = axis_y0 + axis_len_y, linewidth = 2, arrow = axis_arrow) +
  annotate("text", x = axis_x0 + axis_len_x / 2, y = axis_y0 - 0.05 * y_range, label = "UMAP 1", size = 8, family = "Arial", vjust = 1) +
  annotate("text", x = axis_x0 - 0.05 * x_range, y = axis_y0 + axis_len_y / 2, label = "UMAP 2", size = 8, family = "Arial", angle = 90)

p1
ggsave("Figure6a_umap_cluster.pdf", p1, device = cairo_pdf, width = 10, height = 10, units = "in")

crop_colors <- c(Soybean = "#BBDED6", Rice = "#D8BFD8", Wheat = "#FFDAB9")

p2 <- DimPlot(combined_seurat_object, group.by = "crop", reduction = "umap", label = FALSE, pt.size = 1, raster = FALSE, alpha = 1, cols = crop_colors) +
  labs(title = "UMAP by crop") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 40, family = "Arial"),
    plot.title.position = "panel", axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
    legend.text = element_text(size = 20, family = "Arial"), legend.key.height = unit(2, "lines"),
    legend.key.width = unit(1.5, "lines"), legend.position = "right", panel.grid = element_blank(),
    plot.margin = margin(20, 20, 120, 120, "pt")
  ) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 8))) +
  annotate("segment", x = axis_x0, xend = axis_x0 + axis_len_x, y = axis_y0, yend = axis_y0, linewidth = 2, arrow = axis_arrow) +
  annotate("segment", x = axis_x0, xend = axis_x0, y = axis_y0, yend = axis_y0 + axis_len_y, linewidth = 2, arrow = axis_arrow) +
  annotate("text", x = axis_x0 + axis_len_x / 2, y = axis_y0 - 0.05 * y_range, label = "UMAP 1", size = 8, family = "Arial", vjust = 1) +
  annotate("text", x = axis_x0 - 0.05 * x_range, y = axis_y0 + axis_len_y / 2, label = "UMAP 2", size = 8, family = "Arial", angle = 90)

p2
ggsave("Figure6a_umap_crop.pdf", p2, device = cairo_pdf, width = 11, height = 10, units = "in", limitsize = FALSE)

p2_split <- DimPlot(
  combined_seurat_object, group.by = "crop", reduction = "umap", split.by = "crop",
  label = FALSE, pt.size = 2, raster = FALSE, alpha = 1, cols = crop_colors, ncol = 3
) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), panel.grid = element_blank(),
    strip.text = element_text(size = 28, family = "Arial", face = "bold"), strip.background = element_blank(),
    legend.position = "none", plot.margin = margin(20, 20, 120, 120, "pt")
  ) +
  coord_cartesian(xlim = c(x_min - 0.1 * x_range, x_max + 0.02 * x_range), ylim = c(y_min - 0.1 * y_range, y_max + 0.02 * y_range)) +
  annotate("segment", x = axis_x0, xend = axis_x0 + axis_len_x, y = axis_y0, yend = axis_y0, linewidth = 2, arrow = axis_arrow) +
  annotate("segment", x = axis_x0, xend = axis_x0, y = axis_y0, yend = axis_y0 + axis_len_y, linewidth = 2, arrow = axis_arrow) +
  annotate("text", x = axis_x0 + axis_len_x / 2, y = axis_y0 - 0.05 * y_range, label = "UMAP 1", size = 8, family = "Arial", vjust = 1) +
  annotate("text", x = axis_x0 - 0.05 * x_range, y = axis_y0 + axis_len_y / 2, label = "UMAP 2", size = 8, family = "Arial", angle = 90)

p2_split
ggsave("FigureS_umap_crop_split.pdf", p2_split, device = cairo_pdf, width = 24, height = 10, units = "in", limitsize = FALSE)

legend_colors <- c(
  Soybean_1 = "#0FA3B1", Soybean_2 = "#4EC5C1", Soybean_3 = "#B8E3E0",
  Wheat_1 = "#E36414", Wheat_2 = "#F4A261", Wheat_3 = "#FFD6A5",
  Rice_1 = "#6A4C93", Rice_2 = "#9D79BC", Rice_3 = "#D0BDF4"
)

print(table(combined_seurat_object$sample))
combined_seurat_object$sample <- factor(combined_seurat_object$sample, levels = names(legend_colors))

p3 <- DimPlot(
  combined_seurat_object, group.by = "sample", reduction = "umap", label = FALSE,
  pt.size = 1, raster = FALSE, alpha = 1, cols = legend_colors
) +
  labs(title = "UMAP by sample") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 40, family = "Arial"),
    plot.title.position = "panel", axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
    legend.text = element_text(size = 20, family = "Arial"), legend.key.height = unit(2, "lines"),
    legend.key.width = unit(1.5, "lines"), legend.position = "right", panel.grid = element_blank(),
    plot.margin = margin(20, 20, 120, 120, "pt")
  ) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 8))) +
  annotate("segment", x = axis_x0, xend = axis_x0 + axis_len_x, y = axis_y0, yend = axis_y0, linewidth = 2, arrow = axis_arrow) +
  annotate("segment", x = axis_x0, xend = axis_x0, y = axis_y0, yend = axis_y0 + axis_len_y, linewidth = 2, arrow = axis_arrow) +
  annotate("text", x = axis_x0 + axis_len_x / 2, y = axis_y0 - 0.05 * y_range, label = "UMAP 1", size = 8, family = "Arial", vjust = 1) +
  annotate("text", x = axis_x0 - 0.05 * x_range, y = axis_y0 + axis_len_y / 2, label = "UMAP 2", size = 8, family = "Arial", angle = 90)

p3
ggsave("FigureS_umap_sample.pdf", p3, device = cairo_pdf, width = 11, height = 10, units = "in", limitsize = FALSE)

p3_split <- DimPlot(
  combined_seurat_object, group.by = "sample", reduction = "umap", split.by = "sample",
  label = FALSE, pt.size = 3, raster = FALSE, alpha = 1, cols = legend_colors, ncol = 3
) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), panel.grid = element_blank(),
    strip.text = element_text(size = 28, family = "Arial", face = "bold"), strip.background = element_blank(),
    legend.position = "none", plot.margin = margin(20, 20, 120, 120, "pt")
  ) +
  coord_cartesian(xlim = c(x_min - 0.1 * x_range, x_max + 0.02 * x_range), ylim = c(y_min - 0.1 * y_range, y_max + 0.02 * y_range)) +
  annotate("segment", x = axis_x0, xend = axis_x0 + axis_len_x, y = axis_y0, yend = axis_y0, linewidth = 2, arrow = axis_arrow) +
  annotate("segment", x = axis_x0, xend = axis_x0, y = axis_y0, yend = axis_y0 + axis_len_y, linewidth = 2, arrow = axis_arrow) +
  annotate("text", x = axis_x0 + axis_len_x / 2, y = axis_y0 - 0.05 * y_range, label = "UMAP 1", size = 8, family = "Arial", vjust = 1) +
  annotate("text", x = axis_x0 - 0.05 * x_range, y = axis_y0 + axis_len_y / 2, label = "UMAP 2", size = 8, family = "Arial", angle = 90)

p3_split
ggsave("FigureS_umap_sample_split.pdf", p3_split, device = cairo_pdf, width = 25, height = 25, units = "in", limitsize = FALSE)
