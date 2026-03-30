# ===============================================
# Author : ZHENG XINGHAI
# Date   : 2025-09-30
# Project: Soil mscRNA-seq Analysis
# ===============================================

library(Seurat)
library(dplyr)
library(tibble)
library(purrr)

# Get the current working directory and remove extra whitespace
current_dir <- trimws(getwd())
cat("Current dir:", current_dir, "\n")

# Split path and get project root
parts <- strsplit(current_dir, .Platform$file.sep)[[1]]
if (length(parts) >= 2 && all(tail(parts, 2) == c("Soil_Scripts", "Figure4"))) {
  project_root <- paste(parts[1:(length(parts)-2)], collapse = .Platform$file.sep)
} else {
  project_root <- current_dir
}
cat("Project root:", project_root, "\n")

# Define the root directory containing all sample folders
root_dir <- file.path(project_root, "Soil_Matrix")
cat("Root directory for samples:", root_dir, "\n")

# List all sample directories and extract sample names
sample_dirs <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)
sample_names <- basename(sample_dirs)

# Function to process a single sample and compute logCPM expression
process_sample <- function(path) {
  message("Processing: ", path)
  
  adata <- Read10X(path, gene.column = 1)
  seurat_obj <- CreateSeuratObject(counts = adata)
  counts_mat <- seurat_obj@assays$RNA$counts
  
  total_counts <- sum(counts_mat)
  UMI_raw <- rowSums(counts_mat)
  UMI_logCPM <- log(UMI_raw / total_counts * 1e6 + 1)
  
  df <- tibble(
    gene = rownames(seurat_obj),
    expr = UMI_logCPM
  )
  
  df$gene <- gsub("-", "_", df$gene)
  df <- df %>% filter(expr > 0)
  
  return(df)
}

# Apply processing function to all samples and store results in a list
sample_list <- map(sample_dirs, process_sample)
names(sample_list) <- sample_names

# Collect the union of genes across all samples
all_genes <- sample_list %>%
  map(~ .x$gene) %>%
  unlist() %>%
  unique()
message("Total genes in union: ", length(all_genes))

# Initialize expression matrix with genes as rows and samples as columns
expr_mat <- matrix(
  0,
  nrow = length(all_genes),
  ncol = length(sample_list),
  dimnames = list(all_genes, sample_names)
)

# Fill the expression matrix with sample-specific logCPM values
for (i in seq_along(sample_list)) {
  df <- sample_list[[i]]
  expr_mat[df$gene, i] <- df$expr
}

# Export the final logCPM expression matrix as a CSV file
output_file <- "gene_logCPM_matrix.csv"
write.csv(expr_mat, output_file, quote = FALSE)
message("Matrix saved to: ", output_file)