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

current_dir <- trimws(getwd())
setwd(current_dir)
cat("Working directory:", current_dir, "\n")

parts <- strsplit(current_dir, .Platform$file.sep)[[1]]
project_root <- if (length(parts) >= 2 && all(tail(parts, 2) == c("Soil_Scripts", "Figure4"))) paste(parts[1:(length(parts)-2)], collapse = .Platform$file.sep) else current_dir
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

Idents(combined_seurat_object) <- "seurat_clusters"
md <- combined_seurat_object@meta.data
output_df <- data.frame(
  barcode = rownames(md), cluster = as.character(md$seurat_clusters), sample = md$orig.ident,
  crop = md$crop, species = md$species_info, genus = md$genus_info, gene_number = md$gene_number,
  stringsAsFactors = FALSE
)
write.table(output_df, "cell_metadata.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

combined_seurat_object <- PrepSCTFindMarkers(combined_seurat_object)

markers <- FindAllMarkers(combined_seurat_object, assay = "SCT", only.pos = TRUE, logfc.threshold = 0.25)

expression_matrix <- GetAssayData(combined_seurat_object, assay = "SCT", slot = "data")
cluster_ids <- Idents(combined_seurat_object)

markers$cells_expressing <- vapply(seq_len(nrow(markers)), function(i) {
  gene <- markers$gene[i]
  cluster <- markers$cluster[i]
  cells <- names(cluster_ids)[cluster_ids == cluster]
  sum(expression_matrix[gene, cells] > 0)
}, integer(1))

markers <- markers %>% arrange(cluster, desc(avg_log2FC))
write.table(markers, "all_clusters_DEGs.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

counts <- GetAssayData(combined_seurat_object, assay = "SCT", slot = "data")
clusters <- combined_seurat_object$seurat_clusters

markers$gene <- gsub("_", "-", markers$gene)
genes <- unique(markers$gene)
genes <- genes[genes %in% rownames(counts)]

cluster_ids <- sort(unique(clusters))
cluster_names <- paste0("Cluster_", cluster_ids)

fraction_mat <- matrix(0, nrow = length(genes), ncol = length(cluster_ids), dimnames = list(genes, cluster_names))
mean_mat <- fraction_mat

for (cid in cluster_ids) {
  cluster_cells <- names(clusters)[clusters == cid]
  sub_counts <- counts[genes, cluster_cells, drop = FALSE]
  fraction_mat[, paste0("Cluster_", cid)] <- rowSums(sub_counts > 0) / length(cluster_cells)
  mean_mat[, paste0("Cluster_", cid)] <- rowMeans(sub_counts)
}

write.csv(fraction_mat, "gene_fraction_by_cluster.csv", quote = FALSE)
write.csv(mean_mat, "gene_mean_expr_by_cluster.csv", quote = FALSE)

saveRDS(combined_seurat_object, "combined_seurat_object.rds")
combined_seurat_object <- readRDS("combined_seurat_object.rds")

Idents(combined_seurat_object) <- "seurat_clusters"

cluster_annotation <- c(
  "0"="DNA replication|C0", "1"="DNA replication|C1", "2"="Carbohydrate transport|C2",
  "3"="Cell envelope biogenesis|C3", "4"="Carbohydrate metabolism|C4", "5"="Phosphorylation|C5",
  "6"="DNA replication|C6", "7"="Cell envelope biogenesis|C7", "8"="Signal transduction|C8",
  "9"="Antibiotic resistance|C9", "10"="Protein turnover|C10", "11"="Stress response|C11",
  "12"="Protein insertion|C12", "13"="Cell wall remodeling|C13", "14"="Respiration|C14",
  "15"="Protein folding|C15", "16"="DNA repair|C16", "17"="DNA replication|C17",
  "18"="DNA segregation|C18", "19"="Carbohydrate metabolism|C19", "20"="Protein folding|C20",
  "21"="DNA segregation|C21", "22"="Stress response|C22", "23"="DNA segregation|C23",
  "24"="Protein degradation|C24", "25"="Phosphate transport|C25"
)

clusters <- levels(Idents(combined_seurat_object))
cluster_colors <- setNames(colorRampPalette(brewer.pal(12, "Paired"))(length(clusters)), clusters)

umap <- Embeddings(combined_seurat_object, "umap") %>% as.data.frame() %>% rownames_to_column("cell")
umap$cluster <- as.character(Idents(combined_seurat_object)[umap$cell])
umap$cluster_label <- cluster_annotation[umap$cluster]
umap$cluster_label <- factor(umap$cluster_label, levels = cluster_annotation)

set.seed(1024)
umap <- umap[sample(nrow(umap)), ]

x_min <- min(umap$umap_1); x_max <- max(umap$umap_1); y_min <- min(umap$umap_2); y_max <- max(umap$umap_2)
x_range <- x_max - x_min; y_range <- y_max - y_min
axis_len_x <- 0.2 * x_range; axis_len_y <- 0.2 * y_range
axis_x0 <- x_min - 0.05 * x_range; axis_y0 <- y_min - 0.05 * y_range
axis_arrow <- arrow(type = "closed", length = unit(5, "mm"))

axis_layers <- list(
  annotate("segment", x = axis_x0, xend = axis_x0 + axis_len_x, y = axis_y0, yend = axis_y0, linewidth = 2, arrow = axis_arrow),
  annotate("segment", x = axis_x0, xend = axis_x0, y = axis_y0, yend = axis_y0 + axis_len_y, linewidth = 2, arrow = axis_arrow),
  annotate("text", x = axis_x0 + axis_len_x / 2, y = axis_y0 - 0.05 * y_range, label = "UMAP 1", size = 8, family = "Arial", vjust = 1),
  annotate("text", x = axis_x0 - 0.05 * x_range, y = axis_y0 + axis_len_y / 2, label = "UMAP 2", size = 8, family = "Arial", angle = 90, hjust = 0.5)
)

umap_theme <- theme_minimal() + theme(
  plot.title = element_text(hjust = 0.5, size = 40, family = "Arial"),
  plot.title.position = "panel",
  axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
  panel.grid = element_blank(), plot.margin = margin(20, 20, 120, 120, unit = "pt")
)

p1 <- ggplot(umap, aes(umap_1, umap_2, color = cluster_label)) +
  geom_point(size = 0.1, stroke = 0) +
  scale_color_manual(values = setNames(cluster_colors, cluster_annotation)) +
  labs(title = "UMAP by Functional Cluster", color = NULL) +
  umap_theme +
  theme(legend.text = element_text(size = 20, family = "Arial"), legend.key.height = unit(2, "lines")) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 8), title = NULL)) +
  axis_layers

p1
ggsave("Figure4d_cluster.pdf", p1, device = cairo_pdf, width = 14, height = 10, units = "in", limitsize = FALSE)


crop_colors <- c(Soybean = "#BBDED6", Rice = "#D8BFD8", Wheat = "#FFDAB9")

p2 <- DimPlot(
  combined_seurat_object, group.by = "crop", reduction = "umap",
  label = FALSE, pt.size = 0.1, raster = FALSE, alpha = 1, cols = crop_colors
) +
  labs(title = "UMAP by Crop") +
  umap_theme +
  theme(
    legend.text = element_text(size = 20, family = "Arial"),
    legend.key.height = unit(2, "lines"), legend.key.width = unit(1.5, "lines"),
    legend.position = "right"
  ) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 8))) +
  axis_layers

p2
ggsave("FigureS_umap_crop.pdf", p2, device = cairo_pdf, width = 11.5, height = 10, units = "in", limitsize = FALSE)

p2_split <- DimPlot(
  combined_seurat_object, group.by = "crop", reduction = "umap", split.by = "crop",
  label = FALSE, pt.size = 0.1, raster = FALSE, alpha = 1, cols = crop_colors, ncol = 3
) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
    panel.grid = element_blank(), strip.text = element_text(size = 28, family = "Arial", face = "bold"),
    strip.background = element_blank(), legend.position = "none",
    plot.margin = margin(20, 20, 120, 120, unit = "pt")
  ) +
  coord_cartesian(
    xlim = c(x_min - 0.1 * x_range, x_max + 0.02 * x_range),
    ylim = c(y_min - 0.1 * y_range, y_max + 0.02 * y_range)
  ) +
  axis_layers

p2_split
ggsave("FigureS_umap_crop_split.pdf", p2_split, device = cairo_pdf, width = 24, height = 10, units = "in", limitsize = FALSE)


combined_seurat_object$sample_legend <- factor(combined_seurat_object$sample, levels = c(
  "Soybean_1", "Soybean_2", "Soybean_3", "Wheat_1", "Wheat_2", "Wheat_3",
  "Rice_1", "Rice_2", "Rice_3"
))

legend_colors <- c(
  Soybean_1 = "#0FA3B1", Soybean_2 = "#4EC5C1", Soybean_3 = "#B8E3E0",
  Wheat_1 = "#E36414", Wheat_2 = "#F4A261", Wheat_3 = "#FFD6A5",
  Rice_1 = "#6A4C93", Rice_2 = "#9D79BC", Rice_3 = "#D0BDF4"
)

umap$sample_legend <- combined_seurat_object$sample_legend[umap$cell]

p3 <- ggplot(umap, aes(umap_1, umap_2, color = sample_legend)) +
  geom_point(size = 0.1, alpha = 1) +
  scale_color_manual(values = legend_colors) +
  labs(title = "UMAP by Sample", color = NULL) +
  umap_theme +
  theme(
    legend.text = element_text(size = 20, family = "Arial"),
    legend.key.height = unit(2, "lines"), legend.key.width = unit(1.5, "lines"),
    legend.position = "right"
  ) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 8))) +
  axis_layers

p3
ggsave("FigureS_umap_sample.pdf", p3, device = cairo_pdf, width = 11.8, height = 10, units = "in", limitsize = FALSE)

p3_split <- ggplot(umap, aes(umap_1, umap_2, color = sample_legend)) +
  geom_point(size = 0.1, alpha = 1) +
  scale_color_manual(values = legend_colors) +
  facet_wrap(~sample_legend, ncol = 3) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
    panel.grid = element_blank(), strip.text = element_text(size = 24, family = "Arial", face = "bold"),
    strip.background = element_blank(), legend.position = "none",
    plot.margin = margin(20, 20, 100, 100, unit = "pt")
  ) +
  coord_cartesian(
    xlim = c(x_min - 0.1 * x_range, x_max + 0.02 * x_range),
    ylim = c(y_min - 0.1 * y_range, y_max + 0.02 * y_range)
  ) +
  axis_layers

p3_split
ggsave("FigureS_umap_sample_split.pdf", p3_split, device = cairo_pdf, width = 25, height = 25, units = "in", limitsize = FALSE)

Idents(combined_seurat_object) <- "seurat_clusters"

meta <- combined_seurat_object@meta.data

species_info <- vapply(meta$species_info, function(x) {
  if (length(x) == 0 || is.na(x)) NA_character_ else as.character(x)[1]
}, character(1))

genus_info <- vapply(meta$genus_info, function(x) {
  if (length(x) == 0 || is.na(x)) NA_character_ else as.character(x)[1]
}, character(1))

names(species_info) <- rownames(meta)
names(genus_info) <- rownames(meta)

combined_seurat_object$species_info <- species_info
combined_seurat_object$genus_info <- genus_info

species_abundance <- sort(table(species_info, useNA = "no"), decreasing = TRUE)
top_species <- names(species_abundance)[seq_len(min(20, length(species_abundance)))]

combined_seurat_object$top_species <- factor(
  ifelse(combined_seurat_object$species_info %in% top_species, combined_seurat_object$species_info, "Others"),
  levels = c(top_species, "Others")
)

species_colors <- setNames(
  c(colorRampPalette(brewer.pal(12, "Set3"))(length(top_species)), "grey90"),
  levels(combined_seurat_object$top_species)
)

umap_coords <- Embeddings(combined_seurat_object, "umap") %>%
  as.data.frame() %>%
  rownames_to_column("cell")

colnames(umap_coords)[2:3] <- c("UMAP_1", "UMAP_2")
umap_coords$species <- combined_seurat_object$top_species[match(umap_coords$cell, colnames(combined_seurat_object))]
umap_coords <- umap_coords[sample(nrow(umap_coords)), ]

x_min <- min(umap_coords$UMAP_1); x_max <- max(umap_coords$UMAP_1)
y_min <- min(umap_coords$UMAP_2); y_max <- max(umap_coords$UMAP_2)
x_range <- x_max - x_min; y_range <- y_max - y_min
axis_len_x <- 0.2 * x_range; axis_len_y <- 0.2 * y_range
axis_x0 <- x_min - 0.05 * x_range; axis_y0 <- y_min - 0.05 * y_range
axis_arrow <- arrow(type = "closed", length = unit(5, "mm"))

umap_theme <- theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 40, family = "Arial"),
    plot.title.position = "panel",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.text = element_text(size = 15, face = "italic"),
    legend.key.height = unit(1.3, "lines"),
    legend.key.width = unit(1.3, "lines"),
    legend.position = "right",
    panel.grid = element_blank(),
    plot.margin = margin(20, 20, 120, 120, unit = "pt")
  )

p4 <- ggplot(umap_coords, aes(UMAP_1, UMAP_2, color = species)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = species_colors) +
  labs(title = "UMAP by Top 20 Species", color = NULL) +
  umap_theme +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 6))) +
  annotate("segment", x = axis_x0, xend = axis_x0 + axis_len_x, y = axis_y0, yend = axis_y0, linewidth = 2, arrow = axis_arrow) +
  annotate("segment", x = axis_x0, xend = axis_x0, y = axis_y0, yend = axis_y0 + axis_len_y, linewidth = 2, arrow = axis_arrow) +
  annotate("text", x = axis_x0 + axis_len_x / 2, y = axis_y0 - 0.05 * y_range, label = "UMAP 1", size = 8, family = "Arial", vjust = 1) +
  annotate("text", x = axis_x0 - 0.05 * x_range, y = axis_y0 + axis_len_y / 2, label = "UMAP 2", size = 8, family = "Arial", angle = 90, hjust = 0.5)

p4

ggsave(
  "FigureS_umap_species.pdf", p4, device = cairo_pdf,
  width = 13.5, height = 10, units = "in", limitsize = FALSE
)

genus_abundance <- sort(table(genus_info, useNA = "no"), decreasing = TRUE)
top_genus <- names(genus_abundance)[seq_len(min(15, length(genus_abundance)))]

combined_seurat_object$top_genus <- factor(
  ifelse(combined_seurat_object$genus_info %in% top_genus, combined_seurat_object$genus_info, "Others"),
  levels = c(top_genus, "Others")
)

genus_grouped_colors <- setNames(
  c(colorRampPalette(brewer.pal(8, "Dark2"))(length(top_genus)), "grey90"),
  levels(combined_seurat_object$top_genus)
)

umap_coords <- Embeddings(combined_seurat_object, "umap") %>%
  as.data.frame() %>%
  rownames_to_column("cell")

colnames(umap_coords)[2:3] <- c("UMAP_1", "UMAP_2")
umap_coords$genus <- combined_seurat_object$top_genus[match(umap_coords$cell, colnames(combined_seurat_object))]
umap_coords <- umap_coords[sample(nrow(umap_coords)), ]

x_min <- min(umap_coords$UMAP_1); x_max <- max(umap_coords$UMAP_1)
y_min <- min(umap_coords$UMAP_2); y_max <- max(umap_coords$UMAP_2)
x_range <- x_max - x_min; y_range <- y_max - y_min
axis_len_x <- 0.2 * x_range; axis_len_y <- 0.2 * y_range
axis_x0 <- x_min - 0.05 * x_range; axis_y0 <- y_min - 0.05 * y_range

p5 <- ggplot(umap_coords, aes(UMAP_1, UMAP_2, color = genus)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = genus_grouped_colors) +
  labs(title = "UMAP by Top 15 Genus") +
  umap_theme +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 6), title = NULL)) +
  annotate("segment", x = axis_x0, xend = axis_x0 + axis_len_x, y = axis_y0, yend = axis_y0, linewidth = 2, arrow = axis_arrow) +
  annotate("segment", x = axis_x0, xend = axis_x0, y = axis_y0, yend = axis_y0 + axis_len_y, linewidth = 2, arrow = axis_arrow) +
  annotate("text", x = axis_x0 + axis_len_x / 2, y = axis_y0 - 0.05 * y_range, label = "UMAP 1", size = 8, family = "Arial", vjust = 1) +
  annotate("text", x = axis_x0 - 0.05 * x_range, y = axis_y0 + axis_len_y / 2, label = "UMAP 2", size = 8, family = "Arial", angle = 90, hjust = 0.5)

p5

ggsave(
  "FigureS_umap_genus.pdf", p5, device = cairo_pdf,
  width = 11.5, height = 10, units = "in", limitsize = FALSE
)