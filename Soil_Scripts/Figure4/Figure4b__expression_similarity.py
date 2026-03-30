# -*- coding: utf-8 -*-
"""
Created on Thu Aug 21 14:19:46 2025
@author: ZHENG XINGHAI
"""

import pandas as pd
import numpy as np
import itertools
import seaborn as sns
import matplotlib.pyplot as plt

# Configure matplotlib to keep text editable in vector graphics software
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'

# Load the gene-level logCPM expression matrix
df = pd.read_csv("gene_logCPM_matrix.csv", index_col=0)

samples = df.columns.tolist()

# Initialize an empty matrix to store pairwise Pearson correlations between samples
corr_matrix = pd.DataFrame(index=samples, columns=samples, dtype=float)

# Define a function that selects the top 10% highly expressed genes in a given sample
def top10_genes(series):
    threshold = np.percentile(series, 90)
    return series[series >= threshold].index

# Compute pairwise correlations using genes that are highly expressed in both samples
for s1, s2 in itertools.combinations(samples, 2):
    top_s1 = top10_genes(df[s1])
    top_s2 = top10_genes(df[s2])
    inter = top_s1.intersection(top_s2)

    if len(inter) < 3:
        corr = np.nan
    else:
        corr = df.loc[inter, s1].corr(df.loc[inter, s2])

    corr_matrix.loc[s1, s2] = corr
    corr_matrix.loc[s2, s1] = corr

# Set self-correlations to 1 along the diagonal
np.fill_diagonal(corr_matrix.values, 1)

# Generate a heatmap to visualize the sample-to-sample Pearson correlation matrix
plt.figure(figsize=(7, 7))

samples = corr_matrix.index

# Extract species information and simplified sample labels from sample IDs
species = samples.str.split("_").str[0]
simple_labels = samples.str.split("_").str[-1]

ax = sns.heatmap(
    corr_matrix.astype(float),
    cmap="mako",
    vmin=0,
    vmax=1,
    annot=True,
    fmt=".2f",
    square=True,
    linewidths=0.1,
    annot_kws={"size": 15},
    cbar_kws={
        "label": "Pearson Correlation",
        "shrink": 0.8,
        "aspect": 30,
        "pad": 0.02
    }
)

# Adjust the colorbar appearance for better readability
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=15)
cbar.ax.yaxis.label.set_size(18)

# Apply simplified sample labels to the heatmap axes
ax.set_xticks(np.arange(len(simple_labels)) + 0.5)
ax.set_xticklabels(simple_labels, rotation=0, fontsize=15)

ax.set_yticks(np.arange(len(simple_labels)) + 0.5)
ax.set_yticklabels(simple_labels, rotation=0, fontsize=15)

# Add a species annotation bar above the heatmap to highlight crop origins
crop_colors = {
    "Soybean": "#BBDED6",
    "Rice": "#D8BFD8",
    "Wheat": "#FFDAB9"
}

bar_height = 0.04
bar_ax = ax.inset_axes([0, 1.02, 1, bar_height])
bar_ax.axis("off")

num_samples = len(samples)
unique_species = species.unique()

start = 0

for sp in unique_species:

    count = (species == sp).sum()
    width = count / num_samples

    bar_ax.add_patch(
        plt.Rectangle(
            (start / num_samples, 0),
            width,
            1,
            color=crop_colors.get(sp, "grey")
        )
    )

    bar_ax.text(
        start / num_samples + width / 2,
        1.2,
        sp,
        ha="center",
        va="bottom",
        fontsize=16
    )

    start += count

# Save the heatmap as a high-resolution PDF figure
plt.tight_layout()
plt.savefig("Figure4b.pdf", dpi=300, bbox_inches="tight")
plt.show()