# -*- coding: utf-8 -*-
"""
Created on Mon Sep 29 16:50:33 2025
@author: ZHENG XINGHAI
"""

import os
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

anno_file, crop_file, output_dir = "soil_gene_annotation.tsv", "DEGs_by_crop_species_annotated.csv", "COG_enrichment_by_crop"
os.makedirs(output_dir, exist_ok=True)

df_anno = pd.read_csv(anno_file, sep="\t")
df_crop = pd.read_csv(crop_file)
df_crop["gene"] = df_crop["gene"].str.replace("-", "_", regex=False)

cog_dict = {}
for _, row in df_anno.iterrows():
    gene_id, cog_ids = row["query"], str(row["COG_ID"]).strip()
    if cog_ids == "-" or cog_ids.lower() == "nan": continue
    for cog in cog_ids.split(","):
        cog = cog.strip()
        if cog: cog_dict.setdefault(cog, []).append(gene_id)

N = df_anno["query"].nunique()

for crop in sorted(df_crop["cluster"].unique()):
    gene_selected = set(df_crop.loc[df_crop["cluster"] == crop, "gene"].astype(str))
    K = len(gene_selected)
    results = []

    for cog, genes in cog_dict.items():
        genes_set = set(genes)
        n, k = len(genes_set), len(genes_set & gene_selected)
        if k == 0: continue
        results.append([cog, n, k, hypergeom(N, n, K).sf(k - 1)])

    if not results: continue

    df_res = pd.DataFrame(results, columns=["COG_ID", "COG_size", "overlap", "p_value"])
    df_res["q_value"] = multipletests(df_res["p_value"], method="fdr_bh")[1]
    df_res["significant"] = df_res["q_value"] < 0.05
    df_res["rich_factor"] = df_res["overlap"] / df_res["COG_size"]
    df_res = df_res.sort_values("q_value")

    out_file = os.path.join(output_dir, f"COG_enrichment_crop_{crop}.csv")
    df_res.to_csv(out_file, index=False)
    print(f"{crop}: {out_file}")

print(f"All crop enrichment results saved in: {output_dir}")