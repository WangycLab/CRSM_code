# -*- coding: utf-8 -*-
"""
Created on Thu Sep  4 17:39:11 2025
@author: ZHENG XINGHAI
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as mticker

# Configure matplotlib settings to keep text editable in vector graphics software
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'

# Load cell metadata and prepare crop-based color assignments for samples
df_meta = pd.read_csv("cell_metadata.tsv", sep="\t")
crop_colors = {"Soybean": "#BBDED6", "Rice": "#D8BFD8", "Wheat": "#FFDAB9"}
df_meta["crop"] = df_meta["sample"].str.split("_").str[0]
sample_order = sorted(df_meta['sample'].unique())
sample_palette = {sample: crop_colors[sample.split("_")[0]] for sample in sample_order}

# Set seaborn theme for publication-quality plotting
sns.set_theme(style="white", context="talk")

# Calculate the number of cells detected in each sample
cell_counts = df_meta.groupby("sample").size().reset_index(name="cell_count")
cell_counts["sample"] = pd.Categorical(cell_counts["sample"], categories=sample_order, ordered=True)
cell_counts = cell_counts.sort_values("sample")

# Create a two-panel figure containing cell counts and gene counts per cell
fig, (ax_bar, ax_box) = plt.subplots(2, 1, figsize=(7, 9), sharex=True, gridspec_kw={'height_ratios': [1, 1.2]})

# Visualize total cell counts per sample using a bar plot
sns.barplot(x="sample", y="cell_count", data=cell_counts, order=sample_order, palette=sample_palette,
            width=0.55, edgecolor="black", linewidth=1.2, ax=ax_bar)
ax_bar.set_ylabel("Cell Counts", fontsize=20)
ax_bar.set_xlabel("")
ax_bar.tick_params(axis="x", labelbottom=False)
ax_bar.tick_params(axis="y", labelsize=15, width=1.2)
max_cell = cell_counts["cell_count"].max()
upper_bar = (int(max_cell / 5000) + 1) * 5000
ax_bar.set_ylim(0, upper_bar)
ax_bar.yaxis.set_major_locator(mticker.MultipleLocator(5000))
sns.despine(ax=ax_bar)

# Display gene counts per cell using violin plots with boxplot overlays
sns.violinplot(x="sample", y="gene_number", data=df_meta, order=sample_order, palette=sample_palette,
               cut=0, inner=None, linewidth=0, width=0.8, ax=ax_box)
for collection in ax_box.collections:
    collection.set_alpha(0.8)

sns.boxplot(x="sample", y="gene_number", data=df_meta, order=sample_order, width=0.22, showcaps=True,
            boxprops={'facecolor': 'white', 'edgecolor': 'black', 'linewidth': 1.2},
            whiskerprops={'linewidth': 1.2}, capprops={'linewidth': 1.2},
            medianprops={'color': 'black', 'linewidth': 1.8}, showfliers=False, ax=ax_box)

ax_box.set_xlabel("Samples", fontsize=20)
ax_box.set_ylabel("Gene Counts per Cell", fontsize=20)
ax_box.tick_params(axis="x", labelsize=15, rotation=45, width=1.2)
ax_box.tick_params(axis="y", labelsize=15, width=1.2)
ax_box.set_ylim(-10, 200)
ax_box.yaxis.set_major_locator(mticker.MultipleLocator(30))
sns.despine(ax=ax_box)

# Adjust figure layout and export the final figure as a high-resolution PDF
plt.subplots_adjust(hspace=0.15)
plt.tight_layout()
plt.savefig("Figure4a.pdf", dpi=300, bbox_inches='tight')
plt.show()