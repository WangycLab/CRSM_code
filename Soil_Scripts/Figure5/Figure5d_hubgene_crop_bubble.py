# -*- coding: utf-8 -*-
"""
Created on Mon Jan 5 19:13:28 2026
@author: ZHENG XINGHAI
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
import matplotlib

matplotlib.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "font.family": "Arial"
})

fraction_df = pd.read_csv("module_gene_fraction_by_crop.csv", index_col=0)
expr_df = pd.read_csv("module_gene_mean_expr_by_crop.csv", index_col=0)
annot_df = pd.read_csv("combined_module_genes_annotated.csv")

top_genes_df = (annot_df[annot_df["COG_ID"] != "-"]
                .sort_values(["module", "kME"], ascending=[True, False])
                .groupby("module").head(5))
module_genes = top_genes_df["gene_name"].tolist()

print("Number of selected genes:", len(module_genes))
print("Example genes:", module_genes[:5])

fraction_df = fraction_df.loc[fraction_df.index.isin(module_genes)]
expr_df = expr_df.loc[expr_df.index.isin(module_genes)]
if fraction_df.empty:
    raise ValueError("No genes passed the COG_ID and kME filtering criteria.")

scale01 = lambda x: (x - x.min()) / (x.max() - x.min()) if x.max() != x.min() else 0
scaled_fraction_df = fraction_df.apply(scale01, axis=1)
scaled_expr_df = expr_df.apply(scale01, axis=1)

clust_df = scaled_fraction_df.dropna()
gene_order = clust_df.index[leaves_list(linkage(pdist(clust_df.values, metric="euclidean"), method="average"))].tolist()

fraction_melt = scaled_fraction_df.reset_index().melt(id_vars="index", var_name="cluster", value_name="S. Exp. R").rename(columns={"index": "gene"})
expr_melt = scaled_expr_df.reset_index().melt(id_vars="index", var_name="cluster", value_name="S. Exp. L").rename(columns={"index": "gene"})
merged = pd.merge(fraction_melt, expr_melt, on=["gene", "cluster"])

crop_order = ["Rice", "Wheat", "Soybean"]
merged["cluster"] = pd.Categorical(merged["cluster"], categories=crop_order, ordered=True)
merged["gene"] = pd.Categorical(merged["gene"], categories=gene_order, ordered=True)
merged = merged.sort_values(["gene", "cluster"])
merged[["S. Exp. R", "S. Exp. L"]] = merged[["S. Exp. R", "S. Exp. L"]].clip(0, 1)
merged["crop_position"] = merged["cluster"].map({"Rice": 0.8, "Wheat": 1.0, "Soybean": 1.2})

fig, ax = plt.subplots(figsize=(3, 15))
sns.scatterplot(data=merged, x="crop_position", y="gene", size="S. Exp. R", hue="S. Exp. L",
                sizes=(10, 300), palette="flare", edgecolor="silver", linewidth=0, ax=ax)

ax.set(xlabel="Crop", ylabel="Module Gene", xlim=(0.7, 1.3))
ax.set_xticks([0.7, 1.0, 1.3])
ax.set_xticklabels(crop_order, fontsize=15)
ax.tick_params(axis="y", labelsize=15)
for label in ax.get_yticklabels():
    label.set_fontstyle("italic")

ax.legend(title="", fontsize=13, title_fontsize=15, bbox_to_anchor=(2.05, 0), loc="lower right", frameon=False)
plt.subplots_adjust(right=0.75, left=0.15)
plt.savefig("Figure5d.pdf", dpi=300, bbox_inches="tight")
plt.show()