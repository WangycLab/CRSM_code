# -*- coding: utf-8 -*-
"""
Created on Mon Sep 29 16:50:33 2025
@author: ZHENG XINGHAI
"""

import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
import os

# Define input annotation and module gene files and create output folder
anno_file = "soil_gene_annotation.tsv" 
module_file = "combined_module_genes_annotated.csv" 
output_dir = "COG_enrichment_by_module" 
os.makedirs(output_dir, exist_ok=True)

# Load annotation and module gene data
df_anno = pd.read_csv(anno_file, sep="\t")
df_module = pd.read_csv(module_file)

# Standardize module gene names by replacing hyphens with underscores
df_module["gene"] = df_module["gene"].str.replace("-", "_", regex=False)

# Build a dictionary mapping each COG ID to its associated genes
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

# Compute total number of genes in the background for hypergeometric test
N = df_anno["query"].nunique()

# Perform COG enrichment analysis for each co-expression module
modules = sorted(df_module["module"].unique())
for module in modules:

    # Get genes belonging to the current module
    gene_selected = set(df_module[df_module["module"] == module]["gene"].astype(str).tolist())
    K = len(gene_selected)

    results = []
    for cog, genes in cog_dict.items():
        genes_set = set(genes)
        n = len(genes_set) 
        k = len(genes_set & gene_selected) 
        if k == 0:
            continue
        # Calculate hypergeometric p-value for enrichment
        rv = hypergeom(N, n, K)
        pval = rv.sf(k - 1)
        results.append([cog, n, k, pval])

    if results:
        # Compile results into a dataframe
        df_res = pd.DataFrame(results, columns=["COG_ID", "COG_size", "overlap", "p_value"])

        # Correct for multiple testing using FDR (Benjamini-Hochberg)
        reject, pvals_corrected, _, _ = multipletests(df_res["p_value"], method="fdr_bh")
        df_res["q_value"] = pvals_corrected
        df_res["significant"] = reject

        # Compute rich factor (proportion of module genes in each COG)
        df_res["rich_factor"] = df_res["overlap"] / df_res["COG_size"]

        # Sort results by adjusted q-value
        df_res = df_res.sort_values("q_value")

        # Save enrichment results for the current module
        out_file = os.path.join(output_dir, f"COG_enrichment_module_{module}.csv")
        df_res.to_csv(out_file, index=False)
        print(f"module {module} results saved: {out_file}")

# Print completion message after all modules are processed
print(f"All module enrichment results saved in folder: {output_dir}")