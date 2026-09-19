# -*- coding: utf-8 -*-
"""
Created on Mon Jan 5 19:13:28 2026
@author: ZHENG XINGHAI
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, leaves_list

plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "Arial"})

X_LABEL_SIZE, Y_LABEL_SIZE = 30, 30
X_TICK_SIZE, Y_TICK_SIZE = 15, 15

fraction_df = pd.read_csv("gene_fraction_by_crop.csv", index_col=0)
expr_df = pd.read_csv("gene_mean_expr_by_crop.csv", index_col=0)
annot_df = pd.read_csv("DEGs_by_crop_species_annotated.csv")

annot_df = annot_df[~annot_df["COG_category"].isin(["-"]) & ~annot_df["COG_category"].str.contains(",", na=False)]
annot_df["gene_key"] = annot_df["gene"].str.replace("_", "-", regex=False)

common_genes = sorted(set(annot_df["gene_key"]) & set(fraction_df.index) & set(expr_df.index))
fraction_df, expr_df = fraction_df.loc[common_genes], expr_df.loc[common_genes]
annot_df = annot_df[annot_df["gene_key"].isin(common_genes)]

scaled_fraction_df = fraction_df.apply(lambda x: (x - x.min()) / (x.max() - x.min()) if x.max() != x.min() else 0, axis=1)
scaled_expr_df = expr_df.apply(lambda x: (x - x.min()) / (x.max() - x.min()) if x.max() != x.min() else 0, axis=1)

gene_order = scaled_fraction_df.index[leaves_list(linkage(scaled_fraction_df.values, method="average", metric="euclidean"))]

fraction_long = scaled_fraction_df.reset_index().rename(columns={"index": "gene_key"}).melt(id_vars="gene_key", var_name="cluster", value_name="S. Exp. R")
expr_long = scaled_expr_df.reset_index().rename(columns={"index": "gene_key"}).melt(id_vars="gene_key", var_name="cluster", value_name="S. Exp. L")
plot_df = fraction_long.merge(expr_long, on=["gene_key", "cluster"])

top_genes = set()
for clust in plot_df["cluster"].unique():
    top_genes.update(plot_df.loc[plot_df["cluster"] == clust].groupby("gene_key")["S. Exp. R"].max().nlargest(10).index)

plot_df = plot_df[plot_df["gene_key"].isin(top_genes)]
gene_order_final = [g for g in gene_order if g in top_genes]
plot_df["gene_key"] = pd.Categorical(plot_df["gene_key"], categories=gene_order_final, ordered=True)

crop_order = ["Rice", "Wheat", "Soybean"]
plot_df["cluster"] = pd.Categorical(plot_df["cluster"], categories=crop_order, ordered=True)
plot_df["crop_position"] = plot_df["cluster"].map({"Rice": 0.8, "Wheat": 1.0, "Soybean": 1.2})

fig, ax1 = plt.subplots(figsize=(3, 15))

sns.scatterplot(
    data=plot_df, x="crop_position", y="gene_key", size="S. Exp. R", hue="S. Exp. L",
    sizes=(10, 300), palette="flare", edgecolor="silver", ax=ax1
)

ax1.set_xlabel("Crop", fontsize=X_LABEL_SIZE)
ax1.set_ylabel("Crop-Specific DEGs", fontsize=Y_LABEL_SIZE)
ax1.set_xticks([0.7, 1.0, 1.3])
ax1.set_xticklabels(["Rice", "Wheat", "Soybean"])
ax1.tick_params(axis="x", labelsize=X_TICK_SIZE)
ax1.tick_params(axis="y", labelsize=Y_TICK_SIZE)

for label in ax1.get_yticklabels():
    label.set_fontstyle("italic")

ax1.legend(title="", fontsize=13, title_fontsize=15, bbox_to_anchor=(2.05, 0), loc="lower right", frameon=False)

plt.subplots_adjust(right=0.75, left=0.15)
plt.savefig("Figure6c.pdf", dpi=300, bbox_inches="tight")
plt.show()