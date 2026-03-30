# -*- coding: utf-8 -*-
"""
Created on Sun Aug 31 18:21:01 2025
@author: ZHENG XINGHAI
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import pdist
from matplotlib.colors import TwoSlopeNorm
import matplotlib as mpl

# Configure PDF/figure fonts for Adobe Illustrator compatibility
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams['font.family'] = 'Arial'

# Load annotated DEG table
file = "all_clusters_DEGs_annotated.csv"
df = pd.read_csv(file, index_col="gene")

# Extract valid COG IDs for each cluster
cluster_cog_dict = {}
for cluster, group in df.groupby("cluster"):
    cog_list = group["COG_ID"].dropna()
    cog_list = [c for c in cog_list if c != "-"]
    if len(cog_list) > 0:
        cluster_cog_dict[cluster] = cog_list

# Build cluster-by-COG count matrix
all_cogs = sorted(set(cog for lst in cluster_cog_dict.values() for cog in lst))
matrix = pd.DataFrame(0, index=cluster_cog_dict.keys(), columns=all_cogs)
for cluster, cogs in cluster_cog_dict.items():
    for cog in cogs:
        matrix.loc[cluster, cog] += 1

# Compute Pearson correlation between clusters
corr_matrix = matrix.T.corr(method="pearson").fillna(0)

# Perform hierarchical clustering on clusters
linkage_matrix = linkage(pdist(corr_matrix), method="average")
threshold = 1.2
cluster_groups = fcluster(linkage_matrix, t=threshold, criterion="distance")
unique_groups = np.unique(cluster_groups)

# Assign colors to cluster groups for visualization
palette = sns.color_palette("Set3", len(unique_groups))
group_color_map = {grp: palette[i] for i, grp in enumerate(unique_groups)}
row_colors = [group_color_map[grp] for grp in cluster_groups]

# Set Seaborn global plotting style
sns.set_theme(style="ticks", context="talk", font_scale=1.0)

# Define diverging colormap and normalization for heatmap
cmap = sns.diverging_palette(240, 10, as_cmap=True)
norm = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)

# Draw clustered heatmap with dendrograms and row/column colors
g = sns.clustermap(
    corr_matrix,
    row_linkage=linkage_matrix,
    col_linkage=linkage_matrix,
    cmap=cmap,
    norm=norm,
    figsize=(15, 15),
    linewidths=0,
    row_colors=row_colors,
    col_colors=row_colors,
    cbar_kws={"shrink": 0.3},
    annot=False,
    tree_kws={"linewidths": 0.01}
)

# Overlay dendrograms with threshold coloring
dendrogram(linkage_matrix, color_threshold=threshold, ax=g.ax_row_dendrogram, orientation="left", no_labels=True)
g.ax_row_dendrogram.invert_yaxis()
dendrogram(linkage_matrix, color_threshold=threshold, ax=g.ax_col_dendrogram, orientation="top", no_labels=True)

# Adjust heatmap tick labels and colorbar appearance
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), fontsize=25, rotation=0)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=25, rotation=0)
colorbar = g.cax
colorbar.tick_params(labelsize=30)
colorbar.set_title("PCC", fontsize=30, pad=20)

# Adjust thickness and spacing of row/column color annotation bars
g.ax_row_colors.set_position([g.ax_row_colors.get_position().x0 - 0.005, g.ax_row_colors.get_position().y0, g.ax_row_colors.get_position().width * 0.5, g.ax_row_colors.get_position().height])
g.ax_col_colors.set_position([g.ax_col_colors.get_position().x0, g.ax_col_colors.get_position().y0 + 0.015, g.ax_col_colors.get_position().width, g.ax_col_colors.get_position().height * 0.5])

# Save clustered correlation heatmap as PDF
plt.savefig("FigureS4a.pdf", dpi=300, bbox_inches="tight")
plt.show()

# Extract cluster order from hierarchical clustering
reordered_indices = g.dendrogram_row.reordered_ind
sorted_clusters = [matrix.index[i] for i in reordered_indices]
print(sorted_clusters)

# Reload annotated DEGs and rebuild cluster-by-COG matrix
file = "all_clusters_DEGs_annotated.csv"
df = pd.read_csv(file, index_col="gene")
cluster_cog_dict = {}
for cluster, group in df.groupby("cluster"):
    cog_list = group["COG_ID"].dropna()
    cog_list = [c for c in cog_list if c != "-"]
    if len(cog_list) > 0:
        cluster_cog_dict[cluster] = cog_list
all_cogs = sorted(set(cog for lst in cluster_cog_dict.values() for cog in lst))
matrix = pd.DataFrame(0, index=cluster_cog_dict.keys(), columns=all_cogs, dtype=float)
for cluster, cogs in cluster_cog_dict.items():
    for cog in cogs:
        matrix.loc[cluster, cog] += 1

# Remove clusters with zero total counts and recompute correlation and clustering
matrix = matrix.loc[matrix.sum(axis=1) > 0]
corr_matrix = matrix.T.corr(method="pearson")
linkage_matrix = linkage(pdist(corr_matrix), method="average")
threshold = 1.2
cluster_groups = fcluster(linkage_matrix, t=threshold, criterion="distance")
unique_groups = np.unique(cluster_groups)

# Build group-level COG counts and map clusters to groups
group_cog_counts = pd.DataFrame(0, index=unique_groups, columns=all_cogs, dtype=float)
group_cluster_dict = {grp: [] for grp in unique_groups}
for cluster_idx, grp in zip(matrix.index, cluster_groups):
    group_cluster_dict[grp].append(cluster_idx)
    for cog in cluster_cog_dict[cluster_idx]:
        group_cog_counts.loc[grp, cog] += 1

# Compute COG specificity for each group
specificity = pd.DataFrame(index=group_cog_counts.index, columns=group_cog_counts.columns, dtype=float)
for grp in unique_groups:
    in_group = group_cog_counts.loc[grp]
    out_group = group_cog_counts.drop(grp).mean(axis=0)
    specificity.loc[grp] = in_group - out_group

# Assign colors to each group for table visualization
palette = sns.color_palette("hls", len(unique_groups))
group_color_map = {grp: palette[i] for i, grp in enumerate(unique_groups)}

# Extract top 3 most specific COGs per group
top_n = 3
top_cog_per_group_all = specificity.apply(lambda x: x.nlargest(top_n).index.tolist(), axis=1)
top_spec_values_all = specificity.apply(lambda x: x.nlargest(top_n).values.tolist(), axis=1)

# Assemble summary table with group info, colors, clusters, and top COGs
table_df = pd.DataFrame({
    "Group": unique_groups,
    "Color": [group_color_map[grp] for grp in unique_groups],
    "Clusters": [", ".join(str(c) for c in group_cluster_dict[grp]) for grp in unique_groups],
    "Top1_COG_ID": [top_cog_per_group_all.loc[grp][0] for grp in unique_groups],
    "Top1_Specificity": [round(top_spec_values_all.loc[grp][0], 3) for grp in unique_groups],
    "Top2_COG_ID": [top_cog_per_group_all.loc[grp][1] for grp in unique_groups],
    "Top2_Specificity": [round(top_spec_values_all.loc[grp][1], 3) for grp in unique_groups],
    "Top3_COG_ID": [top_cog_per_group_all.loc[grp][2] for grp in unique_groups],
    "Top3_Specificity": [round(top_spec_values_all.loc[grp][2], 3) for grp in unique_groups],
})

# Render the group annotation table as a figure
fig, ax = plt.subplots(figsize=(10, len(unique_groups) * 0.4))
ax.axis("off")
colors = table_df["Color"].tolist()
cell_text = table_df[["Group", "Clusters", "Top1_COG_ID", "Top1_Specificity", "Top2_COG_ID", "Top2_Specificity", "Top3_COG_ID", "Top3_Specificity"]].values
table = ax.table(cellText=cell_text,
                 colLabels=["Group", "Clusters", "Top1 COG_ID", "Top1 Specificity", "Top2 COG_ID", "Top2 Specificity", "Top3 COG_ID", "Top3 Specificity"],
                 colWidths=[0.05, 0.15, 0.12, 0.13, 0.12, 0.13, 0.12, 0.13],
                 cellLoc="center", rowLoc="center", loc="center")

# Apply background colors to group identifier column
for i, color in enumerate(colors):
    table[i + 1, 0].set_facecolor(color)
table.auto_set_font_size(False)
table.set_fontsize(12)
table.scale(1, 1.5)
plt.tight_layout()

# Save group annotation table figure
fig.savefig("FigureS4b.pdf", dpi=300, bbox_inches="tight")
plt.show()

# Export group annotation table to CSV
table_df.to_csv("group_annotation.csv", index=False)