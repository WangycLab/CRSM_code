# -*- coding: utf-8 -*-
"""
Created on Mon Jan 5 19:13:28 2026
@author: ZHENG XINGHAI
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, leaves_list

# Configure matplotlib settings to keep text editable in vector graphics editors
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'

# Load gene expression fraction, mean expression, and annotation datasets
fraction_df = pd.read_csv("gene_fraction_by_crop.csv", index_col=0)
expr_df = pd.read_csv("gene_mean_expr_by_crop.csv", index_col=0)
annot_df = pd.read_csv("DEGs_by_crop_species_annotated.csv")

# Filter annotation records and standardize gene identifiers
annot_df = annot_df[~annot_df['COG_category'].isin(["-"])]
annot_df = annot_df[~annot_df['COG_category'].str.contains(",", na=False)]
annot_df['gene_key'] = annot_df['gene'].str.replace("_", "-", regex=False)

# Identify genes shared across annotation, fraction, and expression datasets
common_genes = (
    set(annot_df['gene_key'])
    & set(fraction_df.index)
    & set(expr_df.index)
)

fraction_df = fraction_df.loc[sorted(common_genes)]
expr_df = expr_df.loc[sorted(common_genes)]
annot_df = annot_df[annot_df['gene_key'].isin(common_genes)]

# Apply row-wise min–max scaling to normalize expression values
scaled_fraction_df = fraction_df.apply(
    lambda x: (x - x.min()) / (x.max() - x.min()) if x.max() != x.min() else 0,
    axis=1
)

scaled_expr_df = expr_df.apply(
    lambda x: (x - x.min()) / (x.max() - x.min()) if x.max() != x.min() else 0,
    axis=1
)

# Perform hierarchical clustering of genes based on scaled expression fractions
Z = linkage(scaled_fraction_df.values, method="average", metric="euclidean")
gene_order = scaled_fraction_df.index[leaves_list(Z)]

# Convert expression matrices into long format suitable for plotting
fraction_long = (
    scaled_fraction_df
    .reset_index()
    .rename(columns={"index": "gene_key"})
    .melt(id_vars="gene_key", var_name="cluster", value_name="S. Exp. R")
)

expr_long = (
    scaled_expr_df
    .reset_index()
    .rename(columns={"index": "gene_key"})
    .melt(id_vars="gene_key", var_name="cluster", value_name="S. Exp. L")
)

plot_df = fraction_long.merge(
    expr_long, on=["gene_key", "cluster"], how="inner"
)

# Select top genes for each crop group based on the highest scaled fraction values
top_genes = []
for clust in plot_df['cluster'].unique():
    sub = plot_df[plot_df['cluster'] == clust]
    top = (
        sub.groupby('gene_key')['S. Exp. R']
        .max()
        .sort_values(ascending=False)
        .head(10)
        .index
    )
    top_genes.extend(top)

top_genes = sorted(set(top_genes))

# Restrict the plotting dataset to the selected genes
plot_df = plot_df[plot_df['gene_key'].isin(top_genes)]

# Apply clustered gene ordering to the visualization dataset
gene_order_final = [g for g in gene_order if g in top_genes]
plot_df['gene_key'] = pd.Categorical(
    plot_df['gene_key'],
    categories=gene_order_final,
    ordered=True
)

# Define crop order and corresponding positions on the x-axis
crop_order = ["Rice", "Wheat", "Soybean"]
crop_position_map = {"Rice": 0.8, "Wheat": 1.0, "Soybean": 1.2}

plot_df['cluster'] = pd.Categorical(
    plot_df['cluster'],
    categories=crop_order,
    ordered=True
)

plot_df['crop_position'] = plot_df['cluster'].map(crop_position_map)

# Create a bubble plot showing scaled gene expression across crop groups
fig, ax1 = plt.subplots(figsize=(3, 15))

sns.scatterplot(
    data=plot_df,
    x="crop_position",
    y="gene_key",
    size="S. Exp. R",
    hue="S. Exp. L",
    sizes=(10, 300),
    palette="flare",
    edgecolor="silver",
    ax=ax1
)

# Configure axis labels and tick formatting
ax1.set_xlabel("Crop", fontsize=25)
ax1.set_ylabel("Crop-Specific DEGs", fontsize=25)

ax1.set_xticks([0.7, 1.0, 1.3])
ax1.set_xticklabels(["Rice", "Wheat", "Soybean"], fontsize=15)

ax1.tick_params(axis='y', labelsize=15)

for label in ax1.get_yticklabels():
    label.set_fontstyle("italic")

# Adjust legend layout for improved readability
leg = ax1.legend(
    title="",
    fontsize=13,
    title_fontsize=15,
    bbox_to_anchor=(2.05, 0),
    loc="lower right",
    frameon=False
)

# Save the figure as a publication-quality PDF
plt.subplots_adjust(right=0.75, left=0.15)
plt.savefig("Figure6c.pdf", dpi=300, bbox_inches="tight")
plt.show()