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

# Ensure exported figures keep editable text
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['font.family'] = 'Arial'

# Load fraction, mean expression, and annotated module gene data
fraction_df = pd.read_csv("module_gene_fraction_by_crop.csv", index_col=0)
expr_df = pd.read_csv("module_gene_mean_expr_by_crop.csv", index_col=0)
annot_df = pd.read_csv("combined_module_genes_annotated.csv")

# Select top module genes by kME and valid COG annotation
annot_df = annot_df[annot_df['COG_ID'] != "-"]
top_genes_df = (
    annot_df
    .sort_values(['module', 'kME'], ascending=[True, False])
    .groupby('module')
    .head(5)
)
module_genes = top_genes_df['gene_name'].tolist()
print("Number of selected genes:", len(module_genes))
print("Example genes:", module_genes[:5])

# Filter fraction and expression matrices for selected module genes
fraction_df = fraction_df.loc[fraction_df.index.isin(module_genes)]
expr_df = expr_df.loc[expr_df.index.isin(module_genes)]
if fraction_df.shape[0] == 0:
    raise ValueError("No genes passed the COG_ID and kME filtering criteria.")

# Perform row-wise 0-1 scaling for fraction and expression values
scaled_fraction_df = fraction_df.apply(
    lambda x: (x - x.min()) / (x.max() - x.min()) if x.max() != x.min() else 0,
    axis=1
)
scaled_expr_df = expr_df.apply(
    lambda x: (x - x.min()) / (x.max() - x.min()) if x.max() != x.min() else 0,
    axis=1
)

# Cluster genes using hierarchical clustering to determine plotting order
clust_df = scaled_fraction_df.dropna()
dist = pdist(clust_df.values, metric="euclidean")
Z = linkage(dist, method="average")
gene_order = clust_df.index[leaves_list(Z)].tolist()

# Melt fraction and expression data for bubble plot visualization
fraction_melt = scaled_fraction_df.reset_index().melt(
    id_vars='index', var_name='cluster', value_name='S. Exp. R'
).rename(columns={'index': 'gene'})
expr_melt = scaled_expr_df.reset_index().melt(
    id_vars='index', var_name='cluster', value_name='S. Exp. L'
).rename(columns={'index': 'gene'})
merged = pd.merge(fraction_melt, expr_melt, on=['gene', 'cluster'])

# Define crop order for consistent plotting
crop_order = ["Rice", "Wheat", "Soybean"]
merged['cluster'] = pd.Categorical(merged['cluster'], categories=crop_order, ordered=True)
merged['gene'] = pd.Categorical(merged['gene'], categories=gene_order, ordered=True)
merged = merged.sort_values(['gene', 'cluster'])

# Clip scaled values to range [0, 1] for visualization
merged['S. Exp. R'] = merged['S. Exp. R'].clip(0, 1)
merged['S. Exp. L'] = merged['S. Exp. L'].clip(0, 1)

# Define x-axis positions for crops to adjust spacing in the plot
crop_position_map = {"Rice": 0.8, "Wheat": 1.00, "Soybean": 1.2}
merged['crop_position'] = merged['cluster'].map(crop_position_map)

# Create bubble plot showing module gene specificity and expression across crops
fig, ax = plt.subplots(figsize=(3, 15))
sns.scatterplot(
    data=merged,
    x="crop_position",
    y="gene",
    size="S. Exp. R",
    hue="S. Exp. L",
    sizes=(10, 300),
    palette="flare",
    edgecolor="silver",
    linewidth=0,
    ax=ax
)

# Configure axis labels and tick parameters
ax.set_xlabel("Crop", fontsize=25)
ax.set_ylabel("Module Gene", fontsize=25)
ax.set_xticks([0.7, 1.0, 1.3])
ax.set_xticklabels(["Rice", "Wheat", "Soybean"], fontsize=15)
ax.set_xlim(0.7, 1.3)
ax.tick_params(axis='y', labelsize=15)
for label in ax.get_yticklabels():
    label.set_fontstyle('italic')

# Customize legend for size and color mappings
leg = ax.legend(
    title="",
    fontsize=13,
    title_fontsize=15,
    bbox_to_anchor=(2.05, 0),
    loc="lower right",
    frameon=False
)

# Adjust layout and save figure to PDF
plt.subplots_adjust(right=0.75, left=0.15)
plt.savefig("Figure5d.pdf", dpi=300, bbox_inches='tight')
plt.show()