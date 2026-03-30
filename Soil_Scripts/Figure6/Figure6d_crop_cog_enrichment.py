# -*- coding: utf-8 -*-
"""
Created on Mon Sep 29 16:50:33 2025
@author: ZHENG XINGHAI
"""

import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
import os

# Define input files and create the output directory for enrichment results
anno_file = "soil_gene_annotation.tsv" 
crop_file = "DEGs_by_crop_species_annotated.csv" 
output_dir = "COG_enrichment_by_crop"
os.makedirs(output_dir, exist_ok=True)

# Load genome annotation data and crop-specific DEG table
df_anno = pd.read_csv(anno_file, sep="\t")
df_crop = pd.read_csv(crop_file)

# Standardize gene identifiers to match annotation format
df_crop["gene"] = df_crop["gene"].str.replace("-", "_", regex=False)

# Build a dictionary mapping each COG ID to its corresponding genes
cog_dict = {}
for _, row in df_anno.iterrows():
    gene_id = row["query"]
    cog_ids = str(row["COG_ID"]).strip()
    if cog_ids == "-" or cog_ids.lower() == "nan":
        continue
    for cog in cog_ids.split(","):
        cog = cog.strip()
        if cog == "":
            continue
        cog_dict.setdefault(cog, []).append(gene_id)

# Determine the total number of unique genes in the background dataset
N = df_anno["query"].nunique()

# Perform COG enrichment analysis separately for each crop group
crops = sorted(df_crop["cluster"].unique())
for crop in crops:

    # Extract the set of genes associated with the current crop
    gene_selected = set(df_crop[df_crop["cluster"] == crop]["gene"].astype(str).tolist())
    K = len(gene_selected) 

    results = []
    for cog, genes in cog_dict.items():
        genes_set = set(genes)
        n = len(genes_set)  
        k = len(genes_set & gene_selected)   
        if k == 0:
            continue
        rv = hypergeom(N, n, K)
        pval = rv.sf(k - 1)
        results.append([cog, n, k, pval])

    if results:
        # Convert enrichment results into a dataframe
        df_res = pd.DataFrame(results, columns=["COG_ID", "COG_size", "overlap", "p_value"])

        # Apply Benjamini–Hochberg correction for multiple testing
        reject, pvals_corrected, _, _ = multipletests(df_res["p_value"], method="fdr_bh")
        df_res["q_value"] = pvals_corrected
        df_res["significant"] = reject

        # Calculate the enrichment rich factor for each COG category
        df_res["rich_factor"] = df_res["overlap"] / df_res["COG_size"]

        # Rank results based on adjusted p-values
        df_res = df_res.sort_values("q_value")

        # Save the enrichment results for the current crop
        out_file = os.path.join(output_dir, f"COG_enrichment_crop_{crop}.csv")
        df_res.to_csv(out_file, index=False)
        print(f"crop {crop} results saved: {out_file}")

# Report completion of the enrichment analysis for all crop groups
print(f"All crop enrichment results saved in folder: {output_dir}")