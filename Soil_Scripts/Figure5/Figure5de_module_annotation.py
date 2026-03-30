# -*- coding: utf-8 -*-
"""
Created on Fri Aug 22 14:43:00 2025
@author: ZHENG XINGHAI
"""

import pandas as pd

# Load the co-expressed module gene table
deg_file = "combined_module_genes.csv"
deg_df = pd.read_csv(deg_file)

# Standardize gene identifiers by replacing hyphens with underscores
deg_df["gene"] = deg_df["gene_name"].astype(str).str.replace("-", "_", regex=False)

# Load the reference genome annotation table
anno_file = "soil_gene_annotation.tsv"
anno_df = pd.read_csv(anno_file, sep="\t")

# Rename annotation ID column to match the module gene identifiers
anno_df = anno_df.rename(columns={"query": "gene"})

# Merge co-expressed module genes with genome annotations
merged_df = pd.merge(deg_df, anno_df, on="gene", how="left")

# Fill missing annotation values with a placeholder
merged_df = merged_df.fillna("-")

# Export the annotated module gene table to a CSV file
output_file = "combined_module_genes_annotated.csv"
merged_df.to_csv(output_file, index=False)

# Print completion message with output file location
print(f"Annotation completed. Output saved to {output_file}")
