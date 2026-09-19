# -*- coding: utf-8 -*-
"""
Created on Fri Aug 22 14:43:00 2025
@author: ZHENG XINGHAI
"""

import pandas as pd

deg_df = pd.read_csv("DEGs_by_crop_species.tsv", sep="\t")
deg_df["gene"] = deg_df["gene"].astype(str).str.replace("-", "_", regex=False)

anno_df = pd.read_csv("soil_gene_annotation.tsv", sep="\t").rename(columns={"query": "gene"})
merged_df = pd.merge(deg_df, anno_df, on="gene", how="left").fillna("-")

output_file = "DEGs_by_crop_species_annotated.csv"
merged_df.to_csv(output_file, index=False)
print(f"Annotation completed. Output saved to {output_file}")