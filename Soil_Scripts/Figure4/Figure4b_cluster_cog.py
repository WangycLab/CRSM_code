# -*- coding: utf-8 -*-
"""
Created on Fri Aug 22 14:43:00 2025
@author: ZHENG XINGHAI
"""

import pandas as pd

deg = pd.read_csv("all_clusters_DEGs.tsv", sep="\t")
deg["gene"] = deg["gene"].astype(str).str.replace("-", "_", regex=False)

anno = pd.read_csv("soil_gene_annotation.tsv", sep="\t").rename(columns={"query": "gene"})
deg = deg.merge(anno, on="gene", how="left").fillna("-")
deg.to_csv("all_clusters_DEGs_annotated.csv", index=False)

deg = deg[deg["COG_ID"] != "-"].copy()
deg["cells_expressing"] = pd.to_numeric(deg["cells_expressing"])

top_deg = (
    deg.sort_values(["cluster", "cells_expressing"], ascending=[True, False])
       .groupby("cluster", group_keys=False)
       .head(2)
)

cog = pd.read_csv("cog-24.def.tab", sep="\t", header=None, usecols=[0, 2],
                  names=["COG_ID", "COG_Name"], engine="python")
top_deg["COG_Name"] = top_deg["COG_ID"].map(dict(zip(cog["COG_ID"], cog["COG_Name"])))

for cluster, sub in top_deg.groupby("cluster"):
    print(f"\n========== Cluster {cluster} ==========")
    for _, row in sub.iterrows():
        print(f"{row['gene']}\t{row['COG_ID']}\t{row['COG_Name']}")

all_selected_genes = top_deg["gene"].drop_duplicates().tolist()
print(f"\n========== All selected genes (unique) ==========\nTotal genes: {len(all_selected_genes)}")
print("\n".join(all_selected_genes))