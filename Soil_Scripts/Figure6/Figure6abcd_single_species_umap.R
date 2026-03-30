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
if (length(parts) >= 2 && all(tail(parts, 2) == c("Soil_Scripts", "Figure6"))) {
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

# Save the merged Seurat object for downstream analysis
saveRDS(combined_seurat_object, file = "combined_seurat_object.rds")

# Load the combined Seurat object from an RDS file
combined_seurat_object <- readRDS(file = "combined_seurat_object.rds")

set.seed(1024)

# Specify the target microbial species for downstream single-cell analysis
target_species <- "Burkholderia cepacia"

# Ensure that species annotation is available in the metadata
if (!"species_info" %in% colnames(combined_seurat_object@meta.data)) {
  stop("species_info column not found in meta.data. Please check upstream annotation steps.")
}

# Identify all cells belonging to the selected species
species_cells <- rownames(combined_seurat_object@meta.data)[
  combined_seurat_object@meta.data$species_info == target_species
]

# Stop the pipeline if no cells are detected for this species
if (length(species_cells) == 0) {
  stop(paste("No cells found for species:", target_species))
}

# Subset the Seurat object to retain only cells from the target species
combined_seurat_object <- subset(combined_seurat_object, cells = species_cells)

cat("Remaining cells after filtering species", target_species, ":", ncol(combined_seurat_object), "\n")

# Perform PCA using previously identified variable genes
combined_seurat_object <- RunPCA(
  combined_seurat_object,
  features = VariableFeatures(combined_seurat_object)
)

# Define a function to automatically detect the optimal PCA elbow point
find_pca_elbow <- function(
    seurat_obj,
    reduction = "pca",
    ndims = 50,
    axis_title_size = 25,
    axis_text_size = 25,
    axis_font = "sans",
    title_size = 30,
    title_font = "sans"
) {
  
  if (!reduction %in% names(seurat_obj@reductions)) {
    stop("Specified reduction not found in the Seurat object.")
  }
  
  # Extract the standard deviation of the top principal components
  pca_sd <- seurat_obj[[reduction]]@stdev[1:ndims]
  
  x <- 1:ndims
  y <- pca_sd
  x_norm <- (x - min(x)) / (max(x) - min(x))
  y_norm <- (y - min(y)) / (max(y) - min(y))
  
  # Estimate the elbow point based on the maximum distance from the diagonal
  distance <- abs(y_norm - (1 - x_norm)) / sqrt(2)
  elbow_point <- which.max(distance)
  
  plot_df <- data.frame(PC = x, SD = y)
  
  # Generate a PCA elbow plot highlighting the suggested PC number
  p <- ggplot(plot_df, aes(x = PC, y = SD)) +
    geom_point(size = 3, color = "#2E86AB") +
    geom_line(color = "#2E86AB") +
    geom_vline(xintercept = elbow_point, linetype = "dashed", color = "red") +
    annotate("text", x = elbow_point + 1, y = max(y),
             label = paste0("Elbow = PC", elbow_point),
             color = "red", size = 8, hjust = 0) +
    labs(
      title = "PCA Elbow Plot",
      x = "Principal Component",
      y = "Standard Deviation"
    ) +
    theme_bw(base_size = 20) +
    theme(
      plot.title = element_text(
        family = title_font, size = title_size, face = "bold", hjust = 0.5
      ),
      axis.title = element_text(
        family = axis_font, size = axis_title_size, face = "bold"
      ),
      axis.text = element_text(
        family = axis_font, size = axis_text_size
      ),
      panel.grid.major = element_line(color = "grey80", linetype = "dotted"),
      panel.grid.minor = element_blank()
    )
  
  print(p)
  
  cat("Suggested number of principal components:", elbow_point, "\n")
  
  return(elbow_point)
}

# Automatically determine the optimal number of principal components
elbow_pc <- find_pca_elbow(combined_seurat_object)

# Perform UMAP embedding and graph-based clustering using the selected PCs
combined_seurat_object <- combined_seurat_object %>%
  RunUMAP(reduction = "pca", dims = 1:elbow_pc) %>%
  FindNeighbors(reduction = "pca", dims = 1:elbow_pc) %>%
  FindClusters(resolution = 0.1)

cat("UMAP and clustering completed. Final object dimension:", dim(combined_seurat_object), "\n")

# Set SCT as the default assay for downstream analysis
DefaultAssay(combined_seurat_object) <- "SCT"

# Use crop identity as the grouping variable
Idents(combined_seurat_object) <- "crop"

# Prepare the SCT model before performing differential expression analysis
combined_seurat_object <- PrepSCTFindMarkers(combined_seurat_object)

# Identify marker genes across crop groups
markers_by_crop <- FindAllMarkers(
  combined_seurat_object,
  assay = "SCT",
  only.pos = TRUE,
  min.pct = 0.05,
  logfc.threshold = 0.25
)

# Export differential expression results for all crop groups
write.table(
  markers_by_crop,
  file = "DEGs_by_crop_species.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Extract the normalized SCT expression matrix
counts <- GetAssayData(
  combined_seurat_object,
  assay = "SCT",
  slot = "data"
)

# Retrieve crop annotation for each cell
crops <- combined_seurat_object$crop
names(crops) <- Cells(combined_seurat_object)

# Collect genes identified in the differential expression analysis
markers_by_crop$gene <- gsub("_", "-", markers_by_crop$gene)
genes <- unique(markers_by_crop$gene)
genes <- genes[genes %in% rownames(counts)]

cat("Number of DEG genes used for crop statistics:", length(genes), "\n")

# Define crop group names
crop_ids <- sort(unique(crops))
crop_names <- paste0(crop_ids)

# Initialize matrices to store expression fraction and mean expression
fraction_mat_crop <- matrix(
  0,
  nrow = length(genes),
  ncol = length(crop_ids),
  dimnames = list(genes, crop_names)
)

mean_mat_crop <- matrix(
  0,
  nrow = length(genes),
  ncol = length(crop_ids),
  dimnames = list(genes, crop_names)
)

# Calculate expression fraction and average expression for each crop group
for (st in crop_ids) {
  crop_cells <- names(crops[crops == st])
  sub_counts <- counts[genes, crop_cells, drop = FALSE]
  
  fraction_mat_crop[, paste0(st)] <-
    rowSums(sub_counts > 0) / length(crop_cells)
  
  mean_mat_crop[, paste0(st)] <-
    rowMeans(sub_counts)
}

# Save gene expression fraction across crop groups
write.csv(
  fraction_mat_crop,
  file = "gene_fraction_by_crop.csv",
  quote = FALSE
)

# Save mean gene expression values across crop groups
write.csv(
  mean_mat_crop,
  file = "gene_mean_expr_by_crop.csv",
  quote = FALSE
)

# Export cluster assignments and metadata for the filtered cells
Idents(combined_seurat_object) <- "seurat_clusters"

# Retrieve metadata table
md_species <- combined_seurat_object@meta.data

if (!"seurat_clusters" %in% colnames(md_species)) {
  stop("seurat_clusters not found. Please ensure FindClusters() has been run.")
}

# Construct a metadata table including cluster, sample, and crop information
output_df_species <- data.frame(
  barcode = rownames(md_species),
  cluster = as.character(md_species$seurat_clusters),
  sample  = md_species$orig.ident,
  crop    = md_species$crop,
  stringsAsFactors = FALSE
)

# Save metadata table for downstream visualization or integration analysis
write.table(
  output_df_species,
  file = "cell_metadata_species.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Define a color palette for clusters based on the number of identified clusters
clusters <- levels(Idents(combined_seurat_object))
n_clusters <- length(clusters)

cluster_colors <- setNames(
  colorRampPalette(brewer.pal(12, "Paired"))(n_clusters),
  clusters
)

# Extract UMAP coordinates and corresponding cluster identities
umap_coords <- Embeddings(combined_seurat_object, "umap") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell")

umap_coords$cluster <- Idents(combined_seurat_object)[umap_coords$cell]

# Determine coordinate ranges for positioning the custom mini-axis
x_min <- min(umap_coords$umap_1)
x_max <- max(umap_coords$umap_1)
y_min <- min(umap_coords$umap_2)
y_max <- max(umap_coords$umap_2)

x_range <- x_max - x_min
y_range <- y_max - y_min

# Define parameters controlling the length and position of the mini-axis
axis_len_x <- 0.2 * x_range
axis_len_y <- 0.2 * y_range

axis_x0 <- x_min - 0.10 * x_range
axis_y0 <- y_min - 0.10 * y_range

axis_arrow <- arrow(
  type = "closed",
  length = unit(6, "mm")
)

# Generate the UMAP visualization colored by cluster identities
p1 <- DimPlot(
  combined_seurat_object,
  reduction = "umap",
  label = FALSE,
  cols = cluster_colors,
  pt.size = 1,
  raster = FALSE
) +
  labs(title = "UMAP by Cluster") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 40, family = "Arial"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.text = element_text(size = 20, family = "Arial"),
    legend.key.height = unit(2, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 90, 90, unit = "pt")
  ) +
  coord_cartesian(clip = "off") +
  guides(color = guide_legend(
    ncol = 1,
    override.aes = list(size = 10)
  )) +
  annotate(
    "segment",
    x = axis_x0,
    xend = axis_x0 + axis_len_x,
    y = axis_y0,
    yend = axis_y0,
    linewidth = 2,
    arrow = axis_arrow
  ) +
  annotate(
    "segment",
    x = axis_x0,
    xend = axis_x0,
    y = axis_y0,
    yend = axis_y0 + axis_len_y,
    linewidth = 2,
    arrow = axis_arrow
  ) +
  annotate(
    "text",
    x = axis_x0 + axis_len_x / 2,
    y = axis_y0 - 0.05 * y_range,
    label = "UMAP 1",
    size = 8,
    family = "Arial",
    vjust = 1
  ) +
  annotate(
    "text",
    x = axis_x0 - 0.05 * x_range,
    y = axis_y0 + axis_len_y / 2,
    label = "UMAP 2",
    size = 8,
    family = "Arial",
    angle = 90,
    hjust = 0.5
  )

p1

# Export the cluster-based UMAP plot as a publication-quality PDF
ggsave(
  filename = "Figure6a_umap_cluster.pdf",
  plot = p1,
  device = cairo_pdf,
  width = 10,
  height = 10,
  units = "in"
)

# Define color scheme representing different crop types
crop_colors <- c(
  "Soybean" = "#BBDED6",
  "Rice"    = "#D8BFD8",
  "Wheat"   = "#FFDAB9"
)

# Extract UMAP coordinates again for positioning the mini-axis
umap_coords <- Embeddings(combined_seurat_object, "umap") %>%
  as.data.frame() %>%
  rownames_to_column("cell")

# Calculate coordinate ranges used to scale the mini-axis
x_min <- min(umap_coords$umap_1)
x_max <- max(umap_coords$umap_1)
y_min <- min(umap_coords$umap_2)
y_max <- max(umap_coords$umap_2)

x_range <- x_max - x_min
y_range <- y_max - y_min

# Define the mini-axis parameters for the crop-colored UMAP plot
axis_len_x <- 0.2 * x_range
axis_len_y <- 0.2 * y_range

axis_x0 <- x_min - 0.05 * x_range
axis_y0 <- y_min - 0.05 * y_range

axis_arrow <- arrow(type = "closed", length = unit(5, "mm"))

# Randomize plotting order of cells to reduce visual bias from overplotting
umap_coords <- umap_coords[sample(nrow(umap_coords)), ]

# Generate the UMAP visualization colored by crop groups
p2 <- DimPlot(
  combined_seurat_object,
  group.by = "crop",
  reduction = "umap",
  label = FALSE,
  pt.size = 1,
  raster = FALSE,
  alpha = 1,
  cols = crop_colors
) +
  labs(title = "UMAP by crop") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 40, family = "Arial"),
    plot.title.position = "panel",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.text = element_text(size = 20, family = "Arial"),
    legend.key.height = unit(2, "lines"),
    legend.key.width = unit(1.5, "lines"),
    legend.position = "right",
    panel.grid = element_blank(),
    plot.margin = margin(20, 20, 120, 120, unit = "pt")
  ) +
  guides(color = guide_legend(
    ncol = 1,
    override.aes = list(size = 8)
  )) +
  annotate(
    "segment",
    x = axis_x0, xend = axis_x0 + axis_len_x,
    y = axis_y0, yend = axis_y0,
    linewidth = 2,
    arrow = axis_arrow
  ) +
  annotate(
    "segment",
    x = axis_x0, xend = axis_x0,
    y = axis_y0, yend = axis_y0 + axis_len_y,
    linewidth = 2,
    arrow = axis_arrow
  ) +
  annotate(
    "text",
    x = axis_x0 + axis_len_x / 2,
    y = axis_y0 - 0.05 * y_range,
    label = "UMAP 1",
    size = 8,
    family = "Arial",
    vjust = 1
  ) +
  annotate(
    "text",
    x = axis_x0 - 0.05 * x_range,
    y = axis_y0 + axis_len_y / 2,
    label = "UMAP 2",
    size = 8,
    family = "Arial",
    angle = 90,
    hjust = 0.5
  )

p2

# Save the crop-based UMAP visualization as a high-resolution PDF
ggsave(
  filename = "Figure6a_umap_crop.pdf",
  plot = p2,
  device = cairo_pdf,
  width = 11,
  height = 10,
  units = "in",
  limitsize = FALSE
)