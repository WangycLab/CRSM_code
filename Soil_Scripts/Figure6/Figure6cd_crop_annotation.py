# -*- coding: utf-8 -*-
"""
Created on Fri Aug 22 14:43:00 2025
@author: ZHENG XINGHAI
"""

import pandas as pd

# Load the differential expression gene table
deg_file = "DEGs_by_crop_species.tsv"
deg_df = pd.read_csv(deg_file, sep="\t")

# Standardize gene identifiers by replacing hyphens with underscores
deg_df["gene"] = deg_df["gene"].astype(str).str.replace("-", "_", regex=False)

# Load the reference genome functional annotation table
anno_file = "soil_gene_annotation.tsv"
anno_df = pd.read_csv(anno_file, sep="\t")

# Rename the annotation ID column to match the gene identifier field
anno_df = anno_df.rename(columns={"query": "gene"})

# Integrate differential expression results with genome annotations
merged_df = pd.merge(deg_df, anno_df, on="gene", how="left")

# Replace missing annotation entries with a placeholder
merged_df = merged_df.fillna("-")

# Export the annotated differential expression results
output_file = "DEGs_by_crop_species_annotated.csv"
merged_df.to_csv(output_file, index=False)

# Print a completion message indicating successful annotation
print(f"Annotation completed. Output saved to {output_file}")