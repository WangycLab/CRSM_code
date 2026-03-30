# -*- coding: utf-8 -*-
"""
Created on Fri Sep 5 14:33:58 2025
@author: ZHENG XINGHAI
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl

# Configure plotting parameters for vector graphics and set default font
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams['pdf.use14corefonts'] = False
plt.rcParams['font.family'] = 'Arial'

# Load species composition data and remove species with zero abundance
df = pd.read_csv("species_abundance_profile.csv", index_col=0)
data = df[(df != 0).any(axis=1)]

# Identify crop-specific sample groups and compute mean abundance for each crop
soybean_cols = [c for c in data.columns if "Soybean" in c]
wheat_cols = [c for c in data.columns if "Wheat" in c]
rice_cols = [c for c in data.columns if "Rice" in c]

soybean_mean = data[soybean_cols].mean(axis=1)
wheat_mean = data[wheat_cols].mean(axis=1)
rice_mean = data[rice_cols]

average_abundance = pd.concat([soybean_mean, wheat_mean, rice_mean], axis=1).mean(axis=1)

# Select top 30 species by average abundance and aggregate remaining species as "Others"
top30_species = average_abundance.nlargest(30).index
top30_data = data.loc[top30_species]
others_data = data.drop(top30_species).sum(axis=0)
top30_data.loc["Others"] = others_data

# Set seaborn style and font parameters for visualization
sns.set_style("white")
sns.set_context("talk")
font_title = {'fontname':'Arial', 'fontsize':20, 'fontweight':'bold'}
font_xy = {'fontname':'Arial', 'fontsize':20}
font_tick = {'fontname':'Arial', 'fontsize':15}

# Generate color palette for stacked bar plot
cmap_choice = 'tab20b'
num_species = top30_data.shape[0]
cmap = cm.get_cmap(cmap_choice, num_species)
colors = [cmap(i) for i in range(num_species)]
colors[-1] = "#cccccc"

# Plot stacked bar chart showing species relative abundance across samples
plt.figure(figsize=(8, 10))
ax = top30_data.T.plot(
    kind='bar',
    stacked=True,
    width=0.6,
    color=colors,
    edgecolor='none'
)

ax.set_ylabel("Relative Abundance", **font_xy)
ax.set_xlabel("Samples", **font_xy)

ax.tick_params(axis='x', labelrotation=45, labelsize=font_tick['fontsize'])
ax.tick_params(axis='y', labelsize=font_tick['fontsize'])

# Customize legend for species with italics style
handles, labels = ax.get_legend_handles_labels()
ax.legend(
    handles,
    labels,
    bbox_to_anchor=(0.45, -0.5),
    loc='upper center',
    fontsize=12,
    ncol=2,
    handlelength=1,
    handleheight=1,
    columnspacing=0.5,
    labelspacing=0.2,
    title="Species",
    title_fontsize=20,
    prop={'style': 'italic'}
)

plt.tight_layout()
plt.savefig("Figure3d.pdf", dpi=300, bbox_inches="tight")
plt.show()

# Reload species data and recompute average abundance for heatmap analysis
df = pd.read_csv("species_abundance_profile.csv", index_col=0)
data = df[(df != 0).any(axis=1)]

soybean_cols = [c for c in data.columns if "Soybean" in c]
wheat_cols = [c for c in data.columns if "Wheat" in c]
rice_cols = [c for c in data.columns if "Rice" in c]

soybean_mean = data[soybean_cols].mean(axis=1)
wheat_mean = data[wheat_cols].mean(axis=1)
rice_mean = data[rice_cols]

average_abundance = pd.concat([soybean_mean, wheat_mean, rice_mean], axis=1).mean(axis=1)

# Select top 30 species and apply log transformation for heatmap
top30_species = average_abundance.nlargest(30).index
top30_data0 = data.loc[top30_species]
log_data = np.log1p(top30_data0)
log_data = log_data[log_data.std(axis=1) != 0]

# Generate clustered heatmap of species abundance
g = sns.clustermap(
    log_data,
    cmap="viridis",
    figsize=(6, 12),
    metric="euclidean",
    method="average",
    z_score=0,
    col_cluster=False,
    xticklabels=True,
    yticklabels=True,
    linewidths=1,
    linecolor='white'
)

# Adjust positions of dendrograms and colorbar
g.ax_col_dendrogram.set_visible(False)
g.ax_heatmap.set_position([0.25, 0.05, 0.65, 0.9])
g.ax_row_dendrogram.set_position([0.15, 0.05, 0.1, 0.9])
g.cax.set_position([0.05, 0.05, 0.02, 0.2])

# Customize heatmap tick labels and font style
g.ax_heatmap.tick_params(axis='x', labelsize=20)
g.ax_heatmap.tick_params(axis='y', labelsize=20)
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=75)
labels = g.ax_heatmap.get_yticklabels()
for label in labels:
    label.set_fontstyle('italic')
g.ax_heatmap.set_yticklabels(labels, rotation=0)

g.ax_heatmap.set_xlabel("Samples", fontsize=40)
g.ax_heatmap.set_ylabel("Species", fontsize=40)
g.cax.set_position([1.6, 0.75, 0.02, 0.15])

plt.savefig("Figure3e.pdf", dpi=300, bbox_inches='tight')
plt.show()