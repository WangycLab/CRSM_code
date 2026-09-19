# -*- coding: utf-8 -*-
"""
Created on Fri Jul 25 16:22:45 2025
@author: ZHENG XINGHAI
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import jensenshannon
from sklearn.manifold import MDS
from scipy.spatial import ConvexHull
import matplotlib.colors as mcolors
import matplotlib as mpl
import seaborn as sns
import glob
import os

# Settings
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["pdf.use14corefonts"] = False
plt.rcParams["font.family"] = "Arial"

crop_colors = {
    "Soybean": "#BBDED6",
    "Rice": "#D8BFD8",
    "Wheat": "#FFDAB9"
}

# Species composition
all_compositions, sample_names = [], []
file_list = glob.glob("soil_species_filted/*.csv")
print(file_list)

for file in file_list:
    df = pd.read_csv(file)
    composition = df.groupby("name")["new_est_reads"].sum()
    composition /= composition.sum()
    sample_name = os.path.basename(file).replace("_sc_taxonomy_filted.csv", "")
    composition.name = sample_name
    all_compositions.append(composition)
    sample_names.append(sample_name)

super_composition_df = pd.concat(all_compositions, axis=1).fillna(0)[sample_names]
super_composition_df.to_csv("species_abundance_profile.csv")

# Jensen-Shannon similarity
data = super_composition_df.T
samples = data.index
n = len(samples)
similarity_matrix = np.zeros((n, n))

for i in range(n):
    for j in range(n):
        similarity_matrix[i, j] = 1 - jensenshannon(data.iloc[i], data.iloc[j])

similarity_df = pd.DataFrame(similarity_matrix, index=samples, columns=samples)

# MDS
mds = MDS(n_components=2, dissimilarity="precomputed", random_state=42)
coords = mds.fit_transform(1 - similarity_df.values)

crops = samples.str.split("_").str[0]
point_colors = crops.map(crop_colors)

plt.figure(figsize=(4, 4))
plt.scatter(coords[:, 0], coords[:, 1], c=point_colors, s=100, edgecolors="none")

for crop in crop_colors:
    mask = crops == crop
    group_coords = coords[mask]
    if len(group_coords) >= 3:
        hull = ConvexHull(group_coords)
        hull_points = group_coords[hull.vertices]
        plt.fill(hull_points[:, 0], hull_points[:, 1], color=mcolors.to_rgba(crop_colors[crop], alpha=0.5), zorder=0)

plt.xlabel("MDS1", fontsize=16)
plt.ylabel("MDS2", fontsize=16)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)

for crop, color in crop_colors.items():
    plt.scatter([], [], c=color, label=crop, s=100)

legend = plt.legend(title="Crop", loc="lower left", frameon=True, fontsize=10, title_fontsize=15)
legend.get_frame().set_edgecolor("black")
legend.get_frame().set_linewidth(1)

plt.savefig("Figure3b.pdf", dpi=300, bbox_inches="tight")
plt.show()

# Heatmap
simple_labels = samples.str[-1]

plt.figure(figsize=(6, 6))
ax = sns.heatmap(
    similarity_df,
    cmap="flare",
    annot=True,
    fmt=".2f",
    square=True,
    linewidths=0.1,
    cbar=True,
    annot_kws={"size": 13},
    cbar_kws={"shrink": 0.8, "pad": 0.02, "aspect": 30}
)

cbar = ax.collections[0].colorbar
cbar.set_label("Jensen-Shannon Similarity", fontsize=16)
cbar.ax.tick_params(labelsize=15, width=2, length=6)

ax.set_xticks(np.arange(len(simple_labels)) + 0.5)
ax.set_xticklabels(simple_labels, rotation=0, fontsize=14)
ax.set_yticks(np.arange(len(simple_labels)) + 0.5)
ax.set_yticklabels(simple_labels, rotation=0, fontsize=14)

# Crop annotation
bar_ax = ax.inset_axes([0, 1.02, 1, 0.035])
bar_ax.axis("off")

num_samples = len(samples)
start = 0

for crop in crops.unique():
    count = (crops == crop).sum()
    width = count / num_samples
    bar_ax.add_patch(plt.Rectangle((start / num_samples, 0), width, 1, color=crop_colors[crop]))
    bar_ax.text(start / num_samples + width / 2, 1.2, crop, ha="center", va="bottom", fontsize=16)
    start += count

plt.tight_layout()
plt.savefig("Figure3c.pdf", dpi=300, bbox_inches="tight")
plt.show()