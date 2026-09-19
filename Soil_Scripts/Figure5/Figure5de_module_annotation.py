# -*- coding: utf-8 -*-
"""
Created on Fri Aug 22 14:43:00 2025
@author: ZHENG XINGHAI
"""

import pandas as pd

deg_df = pd.read_csv("combined_module_genes.csv")
deg_df["gene"] = deg_df["gene_name"].astype(str).str.replace("-", "_", regex=False)

anno_df = pd.read_csv("soil_gene_annotation.tsv", sep="\t").rename(columns={"query": "gene"})
merged_df = deg_df.merge(anno_df, on="gene", how="left").fillna("-")
merged_df.to_csv("combined_module_genes_annotated.csv", index=False)

print("Annotation completed. Output saved to combined_module_genes_annotated.csv")