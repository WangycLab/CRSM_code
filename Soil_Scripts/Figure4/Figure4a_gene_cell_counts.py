# -*- coding: utf-8 -*-
"""
Created on Thu Sep  4 17:39:11 2025
@author: ZHENG XINGHAI
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as mticker

plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "Arial"})

df_meta = pd.read_csv("cell_metadata.tsv", sep="\t")
crop_colors = {"Soybean": "#BBDED6", "Rice": "#D8BFD8", "Wheat": "#FFDAB9"}
df_meta["crop"] = df_meta["sample"].str.split("_").str[0]

sample_order = sorted(df_meta["sample"].unique())
sample_palette = {s: crop_colors[s.split("_")[0]] for s in sample_order}
median_gene_counts = df_meta.groupby("sample")["gene_number"].median().reindex(sample_order)

sns.set_theme(style="white", context="talk")
fig, ax = plt.subplots(figsize=(7, 5.5))

sns.violinplot(
    data=df_meta, x="sample", y="gene_number", order=sample_order, palette=sample_palette,
    cut=0, inner=None, linewidth=0, width=0.8, ax=ax
)
for c in ax.collections: c.set_alpha(0.8)

sns.boxplot(
    data=df_meta, x="sample", y="gene_number", order=sample_order, width=0.22,
    showcaps=True, showfliers=False, ax=ax,
    boxprops={"facecolor": "white", "edgecolor": "black", "linewidth": 1.2},
    whiskerprops={"linewidth": 1.2}, capprops={"linewidth": 1.2},
    medianprops={"color": "black", "linewidth": 1.8}
)

ax.set(xlabel="Samples", ylabel="Gene Counts per Cell", ylim=(-10, 200))
ax.tick_params(axis="x", labelsize=15, rotation=45, width=1.2)
ax.tick_params(axis="y", labelsize=15, width=1.2)
ax.yaxis.set_major_locator(mticker.MultipleLocator(30))

for i, s in enumerate(sample_order):
    ax.text(i, 190, f"{median_gene_counts.loc[s]:.0f}", ha="center", va="bottom",
            fontsize=16, fontfamily="Arial", color="black")

ax.set_xlabel("Samples", fontsize=20)
ax.set_ylabel("Gene Counts per Cell", fontsize=20)
sns.despine(ax=ax)
plt.tight_layout()
plt.savefig("Figure4a.pdf", dpi=300, bbox_inches="tight")
plt.show()