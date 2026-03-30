# -*- coding: utf-8 -*-
"""
Created on Fri Aug 22 14:43:00 2025
@author: ZHENG XINGHAI
"""

import pandas as pd

# Load DEGs for all clusters
deg_file = "all_clusters_DEGs.tsv"
deg_df = pd.read_csv(deg_file, sep="\t")

# Standardize gene names by replacing hyphens with underscores
deg_df["gene"] = deg_df["gene"].astype(str).str.replace("-", "_", regex=False)

# Load gene annotation table
anno_file = "soil_gene_annotation.tsv"
anno_df = pd.read_csv(anno_file, sep="\t")

# Rename column to match gene identifiers in DEGs
anno_df = anno_df.rename(columns={"query": "gene"})

# Merge DEGs with annotations
merged_df = pd.merge(deg_df, anno_df, on="gene", how="left")

# Fill missing annotations with "-"
merged_df = merged_df.fillna("-")

# Save annotated DEGs to CSV
output_file = "all_clusters_DEGs_annotated.csv"
merged_df.to_csv(output_file, index=False)

# Print completion message
print(f"Annotation completed. Output saved to {output_file}")