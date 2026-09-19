# -*- coding: utf-8 -*-
"""
Created on Mon Sep 29 16:50:33 2025
@author: ZHENG XINGHAI
"""

import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
import os

anno_file, module_file, output_dir = "soil_gene_annotation.tsv", "combined_module_genes_annotated.csv", "COG_enrichment_by_module"
os.makedirs(output_dir, exist_ok=True)

df_anno = pd.read_csv(anno_file, sep="\t")
df_module = pd.read_csv(module_file)
df_module["gene"] = df_module["gene"].astype(str).str.replace("-", "_", regex=False)

cog_dict = {}
for _, row in df_anno.iterrows():
    gene_id, cog_ids = row["query"], str(row["COG_ID"]).strip()
    if cog_ids == "-" or cog_ids.lower() == "nan":
        continue
    for cog in map(str.strip, cog_ids.split(",")):
        if cog:
            cog_dict.setdefault(cog, []).append(gene_id)

N = df_anno["query"].nunique()

for module in sorted(df_module["module"].unique()):
    gene_selected = set(df_module.loc[df_module["module"] == module, "gene"].astype(str))
    K = len(gene_selected)
    results = []

    for cog, genes in cog_dict.items():
        genes_set = set(genes)
        n, k = len(genes_set), len(genes_set & gene_selected)
        if k:
            results.append([cog, n, k, hypergeom.sf(k - 1, N, n, K)])

    if results:
        df_res = pd.DataFrame(results, columns=["COG_ID", "COG_size", "overlap", "p_value"])
        df_res["q_value"] = multipletests(df_res["p_value"], method="fdr_bh")[1]
        df_res["significant"] = df_res["q_value"] < 0.05
        df_res["rich_factor"] = df_res["overlap"] / df_res["COG_size"]
        df_res = df_res.sort_values("q_value")
        out_file = os.path.join(output_dir, f"COG_enrichment_module_{module}.csv")
        df_res.to_csv(out_file, index=False)
        print(f"Module {module} results saved: {out_file}")

print(f"All module enrichment results saved in folder: {output_dir}")