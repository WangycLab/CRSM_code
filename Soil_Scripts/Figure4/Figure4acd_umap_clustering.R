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
if (length(parts) >= 2 && all(tail(parts, 2) == c("Soil_Scripts", "Figure4"))) {
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

# Save the merged Seurat object for downstream analysis
saveRDS(combined_seurat_object, file = "combined_seurat_object.rds")

# Load the combined Seurat object from an RDS file
combined_seurat_object <- readRDS(file = "combined_seurat_object.rds")

# UMAP colored by cluster
clusters <- levels(Idents(combined_seurat_object)) 
n_clusters <- length(clusters) 
cluster_colors <- setNames(colorRampPalette(brewer.pal(12, "Paired"))(n_clusters), clusters) 
umap_coords <- Embeddings(combined_seurat_object, "umap") %>% as.data.frame() %>% tibble::rownames_to_column("cell") 
umap_coords$cluster <- Idents(combined_seurat_object)[umap_coords$cell] 

# compute UMAP range for mini-axis
x_min <- min(umap_coords$umap_1) 
x_max <- max(umap_coords$umap_1) 
y_min <- min(umap_coords$umap_2) 
y_max <- max(umap_coords$umap_2) 
x_range <- x_max - x_min 
y_range <- y_max - y_min 
axis_len_x <- 0.2 * x_range 
axis_len_y <- 0.2 * y_range 
axis_x0 <- x_min - 0.05 * x_range 
axis_y0 <- y_min - 0.05 * y_range 
axis_arrow <- arrow(type = "closed", length = unit(5, "mm")) 

# plot UMAP without cluster label
p1 <- DimPlot(combined_seurat_object, reduction = "umap", label = FALSE, cols = cluster_colors, pt.size = 0.1, raster = FALSE) + 
  labs(title = "UMAP by Cluster") + 
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5, size = 40, family = "Arial"),
        plot.title.position = "panel",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.text = element_text(size = 20, family = "Arial"),
        legend.key.height = unit(2, "lines"),
        panel.grid = element_blank(),
        plot.margin = margin(20, 20, 120, 120, unit = "pt")) + 
  guides(color = guide_legend(ncol = 3, override.aes = list(size = 8))) + 
  annotate("segment", x = axis_x0, xend = axis_x0 + axis_len_x, y = axis_y0, yend = axis_y0, linewidth = 2, arrow = axis_arrow) + 
  annotate("segment", x = axis_x0, xend = axis_x0, y = axis_y0, yend = axis_y0 + axis_len_y, linewidth = 2, arrow = axis_arrow) + 
  annotate("text", x = axis_x0 + axis_len_x / 2, y = axis_y0 - 0.05 * y_range, label = "UMAP 1", size = 8, family = "Arial", vjust = 1) + 
  annotate("text", x = axis_x0 - 0.05 * x_range, y = axis_y0 + axis_len_y / 2, label = "UMAP 2", size = 8, family = "Arial", angle = 90, hjust = 0.5) 

p1

# save UMAP cluster plot as PDF
ggsave(filename = "Figure4c_umap_cluster.pdf", plot = p1, device = cairo_pdf, width = 12, height = 10, units = "in", limitsize = FALSE)

# define crop colors
crop_colors <- c("Soybean" = "#BBDED6", "Rice" = "#D8BFD8", "Wheat" = "#FFDAB9")

# extract UMAP coordinates for crop plot
umap_coords <- Embeddings(combined_seurat_object, "umap") %>% as.data.frame() %>% rownames_to_column("cell")

# compute range for mini-axis
x_min <- min(umap_coords$umap_1) 
x_max <- max(umap_coords$umap_1) 
y_min <- min(umap_coords$umap_2) 
y_max <- max(umap_coords$umap_2) 
x_range <- x_max - x_min 
y_range <- y_max - y_min 

# mini-axis parameters
axis_len_x <- 0.2 * x_range 
axis_len_y <- 0.2 * y_range 
axis_x0 <- x_min - 0.05 * x_range 
axis_y0 <- y_min - 0.05 * y_range 
axis_arrow <- arrow(type = "closed", length = unit(5, "mm")) 

# shuffle rows to reduce overplotting
umap_coords <- umap_coords[sample(nrow(umap_coords)), ]

# plot UMAP colored by crop
p2 <- DimPlot(combined_seurat_object, group.by = "crop", reduction = "umap", label = FALSE, pt.size = 0.1, raster = FALSE, alpha = 1, cols = crop_colors) + 
  labs(title = "UMAP by Crop") + 
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5, size = 40, family = "Arial"),
        plot.title.position = "panel",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.text = element_text(size = 20, family = "Arial"),
        legend.key.height = unit(2, "lines"),
        legend.key.width = unit(1.5, "lines"),
        legend.position = "right",
        panel.grid = element_blank(),
        plot.margin = margin(20, 20, 120, 120, unit = "pt")) + 
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 8))) + 
  annotate("segment", x = axis_x0, xend = axis_x0 + axis_len_x, y = axis_y0, yend = axis_y0, linewidth = 2, arrow = axis_arrow) + 
  annotate("segment", x = axis_x0, xend = axis_x0, y = axis_y0, yend = axis_y0 + axis_len_y, linewidth = 2, arrow = axis_arrow) + 
  annotate("text", x = axis_x0 + axis_len_x / 2, y = axis_y0 - 0.05 * y_range, label = "UMAP 1", size = 8, family = "Arial", vjust = 1) + 
  annotate("text", x = axis_x0 - 0.05 * x_range, y = axis_y0 + axis_len_y / 2, label = "UMAP 2", size = 8, family = "Arial", angle = 90, hjust = 0.5)

p2

# save UMAP crop plot as PDF
ggsave(filename = "FigureS3_umap_crop.pdf", plot = p2, device = cairo_pdf, width = 11.5, height = 10, units = "in", limitsize = FALSE)

# create legend variable for sample
combined_seurat_object$sample_legend <- combined_seurat_object$sample

# define colors for sample legend
legend_colors <- c("Soybean_1" = "#0FA3B1", "Soybean_2" = "#4EC5C1", "Soybean_3" = "#B8E3E0",
                   "Wheat_1" = "#E36414", "Wheat_2" = "#F4A261", "Wheat_3" = "#FFD6A5",
                   "Rice_1" = "#6A4C93", "Rice_2" = "#9D79BC", "Rice_3" = "#D0BDF4")

# ensure factor levels match legend colors
combined_seurat_object$sample_legend <- factor(combined_seurat_object$sample_legend, levels = names(legend_colors))

# extract UMAP coordinates and add sample legend
umap_coords <- Embeddings(combined_seurat_object, "umap") %>% as.data.frame()
colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
umap_coords <- tibble::rownames_to_column(umap_coords, "cell")
umap_coords$sample_legend <- combined_seurat_object$sample_legend[umap_coords$cell]

# shuffle rows to reduce overplotting
umap_coords <- umap_coords[sample(nrow(umap_coords)), ]

# compute UMAP range for mini-axis
x_min <- min(umap_coords$UMAP_1)
x_max <- max(umap_coords$UMAP_1)
y_min <- min(umap_coords$UMAP_2)
y_max <- max(umap_coords$UMAP_2)
x_range <- x_max - x_min
y_range <- y_max - y_min
axis_len_x <- 0.2 * x_range
axis_len_y <- 0.2 * y_range
axis_x0 <- x_min - 0.05 * x_range
axis_y0 <- y_min - 0.05 * y_range
axis_arrow <- arrow(type = "closed", length = unit(5, "mm"))

# plot UMAP colored by sample
p3 <- ggplot(umap_coords, aes(x = UMAP_1, y = UMAP_2, color = sample_legend)) +
  geom_point(size = 0.1, alpha = 1) +
  scale_color_manual(values = legend_colors) +
  labs(title = "UMAP by Sample", color = NULL) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 40, family = "Arial"),
        plot.title.position = "panel",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.text = element_text(size = 20, family = "Arial"),
        legend.key.height = unit(2, "lines"),
        legend.key.width = unit(1.5, "lines"),
        legend.position = "right",
        panel.grid = element_blank(),
        plot.margin = margin(20, 20, 120, 120, unit = "pt")) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 8))) +
  annotate("segment", x = axis_x0, xend = axis_x0 + axis_len_x, y = axis_y0, yend = axis_y0, linewidth = 2, arrow = axis_arrow) +
  annotate("segment", x = axis_x0, xend = axis_x0, y = axis_y0, yend = axis_y0 + axis_len_y, linewidth = 2, arrow = axis_arrow) +
  annotate("text", x = axis_x0 + axis_len_x / 2, y = axis_y0 - 0.05 * y_range, label = "UMAP 1", size = 8, family = "Arial", vjust = 1) +
  annotate("text", x = axis_x0 - 0.05 * x_range, y = axis_y0 + axis_len_y / 2, label = "UMAP 2", size = 8, family = "Arial", angle = 90, hjust = 0.5)

p3

# save UMAP sample plot as PDF
ggsave(filename = "FigureS3_umap_sample.pdf", plot = p3, device = cairo_pdf, width = 11.8, height = 10, units = "in", limitsize = FALSE)

# compute species abundance and select top 20
species_abundance <- combined_seurat_object@meta.data %>%
  group_by(species_info) %>%
  summarise(cell_count = n()) %>%
  arrange(desc(cell_count))
top_species <- head(species_abundance$species_info, 20)

# create top_species column, others are labeled "Others"
combined_seurat_object$top_species <- ifelse(
  combined_seurat_object$species_info %in% top_species,
  combined_seurat_object$species_info,
  "Others"
)
species_levels <- c(top_species, "Others")
combined_seurat_object$top_species <- factor(combined_seurat_object$top_species, levels = species_levels)

# generate colors for top 20 species
top_colors <- colorRampPalette(brewer.pal(12, "Set3"))(20)
species_colors <- setNames(c(top_colors, "grey90"), species_levels)

# extract UMAP coordinates and add species info
umap_coords <- Embeddings(combined_seurat_object, "umap") %>% as.data.frame() %>% tibble::rownames_to_column("cell")
umap_coords$species <- combined_seurat_object$top_species[umap_coords$cell]

# shuffle rows to reduce overplotting
umap_coords <- umap_coords[sample(nrow(umap_coords)), ]

# compute UMAP range for mini-axis
x_min <- min(umap_coords$umap_1)
x_max <- max(umap_coords$umap_1)
y_min <- min(umap_coords$umap_2)
y_max <- max(umap_coords$umap_2)
x_range <- x_max - x_min
y_range <- y_max - y_min
axis_len_x <- 0.2 * x_range
axis_len_y <- 0.2 * y_range
axis_x0 <- x_min - 0.05 * x_range
axis_y0 <- y_min - 0.05 * y_range
axis_arrow <- arrow(type = "closed", length = unit(5, "mm"))

# plot UMAP colored by top 20 species
p4 <- ggplot(umap_coords, aes(x = umap_1, y = umap_2, color = species)) +
  geom_point(size = 0.1, alpha = 1) +
  scale_color_manual(values = species_colors) +
  labs(title = "UMAP by Top 20 Species", color = NULL) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 40, family = "Arial"),
        plot.title.position = "panel",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.text = element_text(size = 15, face = "italic"),
        legend.key.height = unit(1.3, "lines"),
        legend.key.width = unit(1.3, "lines"),
        legend.position = "right",
        panel.grid = element_blank(),
        plot.margin = margin(20, 20, 120, 120, unit = "pt")) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 6))) +
  annotate("segment", x = axis_x0, xend = axis_x0 + axis_len_x, y = axis_y0, yend = axis_y0, linewidth = 2, arrow = axis_arrow) +
  annotate("segment", x = axis_x0, xend = axis_x0, y = axis_y0, yend = axis_y0 + axis_len_y, linewidth = 2, arrow = axis_arrow) +
  annotate("text", x = axis_x0 + axis_len_x / 2, y = axis_y0 - 0.05 * y_range, label = "UMAP 1", size = 8, family = "Arial", vjust = 1) +
  annotate("text", x = axis_x0 - 0.05 * x_range, y = axis_y0 + axis_len_y / 2, label = "UMAP 2", size = 8, family = "Arial", angle = 90, hjust = 0.5)

p4

# save UMAP top 20 species as PDF
ggsave(filename = "FigureS3_umap_species.pdf", plot = p4, device = cairo_pdf, width = 13.5, height = 10, units = "in", limitsize = FALSE)

# loop to highlight each of the top 20 species individually
output_dir <- "UMAP_top20_species"
if (!dir.exists(output_dir)) dir.create(output_dir)

for (sp in top_species) {
  umap_df <- Embeddings(combined_seurat_object, "umap") %>% as.data.frame()
  colnames(umap_df)[1:2] <- c("UMAP1", "UMAP2")
  umap_df$species_info <- combined_seurat_object$species_info
  umap_df$highlight <- ifelse(umap_df$species_info == sp, sp, "Others")
  umap_df <- umap_df[sample(nrow(umap_df)), ]
  
  plot_colors <- c()
  plot_colors[sp] <- species_colors[sp]
  plot_colors["Others"] <- "grey90"
  
  p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
    geom_point(data = subset(umap_df, highlight == "Others"), aes(color = highlight), size = 0.1, alpha = 0.3) +
    geom_point(data = subset(umap_df, highlight == sp), aes(color = highlight), size = 0.1, alpha = 1) +
    scale_color_manual(values = plot_colors) +
    ggtitle(sp) +
    theme_minimal() +
    theme(plot.title = element_text(face = "italic", size = 20, hjust = 0.5),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.grid = element_blank(),
          panel.border = element_blank())
  
  ggsave(file.path(output_dir, paste0("UMAP_", sp, ".jpg")), plot = p, width = 5, height = 5, dpi = 300)
}

# compute genus abundance and select top 15
genus_abundance <- combined_seurat_object@meta.data %>%
  group_by(genus_info) %>%
  summarise(cell_count = n()) %>%
  arrange(desc(cell_count))
top_genus <- head(genus_abundance$genus_info, 15)

# create top15_genus column, others labeled "Others"
combined_seurat_object$top_genus <- ifelse(
  combined_seurat_object$genus_info %in% top_genus,
  combined_seurat_object$genus_info,
  "Others"
)
combined_seurat_object$top_genus <- factor(combined_seurat_object$top_genus, levels = c(top_genus, "Others"))

# generate colors for top 15 genus
top_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(15)
genus_grouped_colors <- setNames(c(top_colors, "grey90"), levels(combined_seurat_object$top_genus))

# extract UMAP coordinates and add genus info
umap_coords <- Embeddings(combined_seurat_object, "umap") %>% as.data.frame() %>% tibble::rownames_to_column("cell")
colnames(umap_coords)[2:3] <- c("UMAP_1", "UMAP_2")
umap_coords$genus <- combined_seurat_object$top_genus[umap_coords$cell]
umap_coords <- umap_coords[sample(nrow(umap_coords)), ]

# compute UMAP range for mini-axis
x_min <- min(umap_coords$UMAP_1)
x_max <- max(umap_coords$UMAP_1)
y_min <- min(umap_coords$UMAP_2)
y_max <- max(umap_coords$UMAP_2)
x_range <- x_max - x_min
y_range <- y_max - y_min
axis_len_x <- 0.2 * x_range
axis_len_y <- 0.2 * y_range
axis_x0 <- x_min - 0.05 * x_range
axis_y0 <- y_min - 0.05 * y_range
axis_arrow <- arrow(type = "closed", length = unit(5, "mm"))

# plot UMAP by top 15 genus
p5 <- ggplot(umap_coords, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = genus), size = 0.1, alpha = 1) +
  scale_color_manual(values = genus_grouped_colors) +
  labs(title = "UMAP by Top 15 Genus") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 40, family = "Arial"),
        plot.title.position = "panel",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.text = element_text(size = 15, face = "italic"),
        legend.key.height = unit(1.5, "lines"),
        panel.grid = element_blank(),
        plot.margin = margin(20, 20, 120, 120, unit = "pt")) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 6), title = NULL)) +
  annotate("segment", x = axis_x0, xend = axis_x0 + axis_len_x, y = axis_y0, yend = axis_y0, linewidth = 2, arrow = axis_arrow) +
  annotate("segment", x = axis_x0, xend = axis_x0, y = axis_y0, yend = axis_y0 + axis_len_y, linewidth = 2, arrow = axis_arrow) +
  annotate("text", x = axis_x0 + axis_len_x / 2, y = axis_y0 - 0.05 * y_range, label = "UMAP 1", size = 8, family = "Arial", vjust = 1) +
  annotate("text", x = axis_x0 - 0.05 * x_range, y = axis_y0 + axis_len_y / 2, label = "UMAP 2", size = 8, family = "Arial", angle = 90, hjust = 0.5)

p5

# save UMAP top 15 genus as PDF
ggsave("FigureS3_umap_genus.pdf", plot = p5, device = cairo_pdf, width = 11.5, height = 10, units = "in", limitsize = FALSE)

# loop to highlight each top 15 genus individually
output_dir <- "UMAP_top15_genus"
if (!dir.exists(output_dir)) dir.create(output_dir)

for (gn in top_genus) {
  umap_df <- Embeddings(combined_seurat_object, "umap") %>% as.data.frame()
  colnames(umap_df)[1:2] <- c("UMAP1", "UMAP2")
  umap_df$genus_info <- combined_seurat_object$genus_info
  umap_df$highlight <- ifelse(umap_df$genus_info == gn, gn, "Others")
  
  plot_colors <- c()
  plot_colors[gn] <- genus_grouped_colors[gn]
  plot_colors["Others"] <- "grey90"
  
  p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
    geom_point(data = subset(umap_df, highlight == "Others"), aes(color = highlight), size = 0.1, alpha = 0.3) +
    geom_point(data = subset(umap_df, highlight == gn), aes(color = highlight), size = 0.1, alpha = 1) +
    scale_color_manual(values = plot_colors) +
    ggtitle(gn) +
    theme_minimal() +
    theme(plot.title = element_text(face = "italic", size = 20, hjust = 0.5),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.grid = element_blank(),
          panel.border = element_blank())
  
  ggsave(file.path(output_dir, paste0("UMAP_", gn, ".jpg")), plot = p, width = 5, height = 5, dpi = 300)
}