# -*- coding: utf-8 -*-
"""
Created on Sun Aug 24 18:32:42 2025
@author: ZHENG XINGHAI
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors

plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none", "font.family": "Arial"})

fraction_df = pd.read_csv("gene_fraction_by_cluster.csv", index_col=0)
expr_df = pd.read_csv("gene_mean_expr_by_cluster.csv", index_col=0)
annot_df = pd.read_csv("all_clusters_DEGs_annotated.csv")
annot_df = annot_df[~annot_df["COG_ID"].isin(["-"]) & ~annot_df["COG_ID"].str.contains(",", na=False)].copy()
annot_df["gene_key"] = annot_df["gene"].str.replace("_", "-", regex=False)

anno = annot_df[["gene_key", "COG_ID"]]
fraction_df = fraction_df.merge(anno, left_index=True, right_on="gene_key", how="inner")
expr_df = expr_df.merge(anno, left_index=True, right_on="gene_key", how="inner")

fraction_melt = fraction_df.melt(id_vars=["gene_key", "COG_ID"], var_name="cluster", value_name="fraction")
expr_melt = expr_df.melt(id_vars=["gene_key", "COG_ID"], var_name="cluster", value_name="expression")
merged = fraction_melt.merge(expr_melt, on=["gene_key", "COG_ID", "cluster"])

selected_genes = [
    "AOX55-RS00005", "AOX55-RS18165", "OYW20-RS00005", "CLDAP-RS00005", "N7L95-RS00365",
    "PX653-RS00005", "OJF58-RS00005", "FDP22-RS00005", "QTJ18-RS01415", "FE840-RS00005",
    "CPT03-RS00005", "G7074-RS00005", "DVA44-RS00005", "GT391-RS00005", "KW115-RS00005",
    "METAL-RS00005", "O4O04-RS17235", "LOC62-05G007437", "D4L85-RS00005", "DTQ70-RS00005",
    "IT6-RS00005", "PG1C-RS00005", "CUN60-RS13030", "C4375-RS00005", "L1P08-RS00005",
    "PY308-RS00005", "BSY16-RS00005", "K9M53-RS00005", "LK994-RS00005", "LUA81-RS00005",
    "MXMO3-RS00005", "CFE-RS00005", "E3U44-RS19160", "J3U87-RS00005", "CCALI-RS00005",
    "D5261-RS00005", "AOP6-RS00040", "Q0X14-RS00005", "BSK21-RS00005", "K0B96-RS00005",
    "JL105-RS00005", "B5M13-RS00005", "EXU85-RS00005", "M6B22-RS00005", "M6D93-RS00005"
]

plot_df = merged[merged["gene_key"].isin(selected_genes)].copy()
gene_cog_map = plot_df[["gene_key", "COG_ID"]].drop_duplicates()
cog_counts = gene_cog_map["COG_ID"].value_counts()
cog_index, unique_map = {}, {}

for _, row in gene_cog_map.iterrows():
    gene, cog = row["gene_key"], row["COG_ID"]
    if cog_counts[cog] > 1: cog_index[cog] = cog_index.get(cog, 0) + 1
    unique_map[gene] = f"{cog} ({cog_index[cog]})" if cog_counts[cog] > 1 else cog

plot_df["COG_label"] = plot_df["gene_key"].map(unique_map)
cog_order = [unique_map[g] for g in selected_genes if g in unique_map]
plot_df["COG_label"] = pd.Categorical(plot_df["COG_label"], categories=cog_order, ordered=True)
plot_df[["fraction", "expression"]] = plot_df[["fraction", "expression"]].clip(upper=1)
plot_df = plot_df.sort_values(["cluster", "COG_label"])

df_meta = pd.read_csv("cell_metadata.tsv", sep="\t")
df_meta["cluster"] = df_meta["cluster"].map(lambda x: f"Cluster_{x}")
count_table = pd.crosstab(df_meta["cluster"], df_meta["genus"])
prop_table = count_table.div(count_table.sum(axis=1), axis=0)

def keep_top_n(row, n=1):
    top = row.nlargest(n).index
    row.loc[~row.index.isin(top)] = 0
    return row

prop_table = prop_table.apply(keep_top_n, axis=1)
prop_table = prop_table.loc[:, (prop_table != 0).any(axis=0)]

cluster_order = [f"Cluster_{i}" for i in range(26)]
prop_table = prop_table.reindex(cluster_order)
genus_list = prop_table.columns.tolist()
palette = sns.color_palette("twilight", n_colors=len(genus_list))
genus_color_map = dict(zip(genus_list, map(mcolors.to_hex, palette)))
print(genus_color_map)

plot_df[["fraction", "expression"]] = plot_df.groupby("cluster")[["fraction", "expression"]].transform(
    lambda x: (x - x.min()) / (x.max() - x.min()) if x.max() != x.min() else 0
)

cluster_annotation = {str(i): f"C{i}" for i in range(26)}
plot_df["cluster"] = plot_df["cluster"].astype(str).str.replace("Cluster_", "", regex=False)
plot_df["cluster_label"] = pd.Categorical(
    plot_df["cluster"].map(cluster_annotation),
    categories=[f"C{i}" for i in range(26)], ordered=True
)

prop_table.index = prop_table.index.astype(str).str.replace("Cluster_", "", regex=False)
prop_table = prop_table.reindex([str(i) for i in range(26)])
prop_table.index = [f"C{i}" for i in range(26)]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(23, 16), sharey=True, gridspec_kw={"width_ratios": [3, 1]})

sns.scatterplot(
    data=plot_df, x="COG_label", y="cluster_label", size="fraction", hue="expression",
    sizes=(10, 300), palette="crest", edgecolor="black", ax=ax1
)

ax1.set(xlabel="COG IDs", ylabel="Functional clusters")
ax1.tick_params(axis="x", rotation=90, labelsize=20)
ax1.tick_params(axis="y", labelsize=20)
ax1.xaxis.label.set_size(30)
ax1.yaxis.label.set_size(30)

handles, labels = ax1.get_legend_handles_labels()
new_handles, new_labels = [], []

for h, l in zip(handles, labels):
    if l == "fraction": new_handles.append(h); new_labels.append("S. Exp. R")
    elif l == "expression": new_handles.append(h); new_labels.append("S. Exp. L")
    elif l not in ["size", "hue"]: new_handles.append(h); new_labels.append(l)

leg = ax1.legend(new_handles, new_labels, title="", ncol=2, fontsize=18, bbox_to_anchor=(0.05, -0.4), loc="upper left")
for text in leg.get_texts(): text.set_fontsize(25 if text.get_text() in ["S. Exp. L", "S. Exp. R"] else 18)

bottom = np.zeros(len(prop_table))
handles, labels = [], []

for genus in prop_table.columns:
    values = prop_table[genus].values
    bar = ax2.barh(prop_table.index, values, left=bottom, color=genus_color_map.get(genus, "#CCCCCC"), height=0.6, edgecolor="none")
    bottom += values
    handles.append(bar[0]); labels.append(genus)

ax2.set_yticks(range(len(prop_table)))
ax2.set_yticklabels(prop_table.index, fontsize=20)
ax2.set_xlabel("Genus ratio", fontsize=30)
ax2.set_xticks([0, 0.25, 0.5, 0.75, 1])
ax2.set_xticklabels(["0", "0.25", "0.5", "0.75", "1.0"], fontsize=20)
ax2.set_ylabel("")
ax2.spines[["top", "right", "bottom"]].set_visible(False)

legend2 = ax2.legend(
    handles, labels, title="Core genus", fontsize=18, title_fontsize=25,
    labelspacing=0.3, ncol=2, columnspacing=1.5, handletextpad=0.5,
    bbox_to_anchor=(-0.1, -0.15), loc="upper left"
)
for text in legend2.get_texts(): text.set_fontstyle("italic")
for patch in legend2.get_patches(): patch.set_edgecolor("none")

plt.tight_layout()
plt.savefig("Figure4d.pdf", dpi=300, bbox_inches="tight")
plt.show()