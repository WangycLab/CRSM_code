# -*- coding: utf-8 -*-
"""
Created on Fri Sep 5 14:33:58 2025
@author: ZHENG XINGHAI
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

# Set basic plotting parameters for all figures
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
point_size = 40
alpha = 0.8
figsize = (3, 3.5)
sns.set(style="ticks", context="paper")

# Load metatranscriptomic and metagenomic species abundance tables
msc = pd.read_csv("msc_species.tsv", sep="\t")
meta = pd.read_csv("meta_speices.tsv", sep="\t")

# Prepare species abundance data
msc_species = (
    msc[["name", "fraction_total_reads"]]
    .set_index("name")
    .rename(columns={"fraction_total_reads": "MscRNA"})
)

meta_species = (
    meta[["name", "fraction_total_reads"]]
    .set_index("name")
    .rename(columns={"fraction_total_reads": "Metagenomic"})
)

# Merge the two datasets on shared species
df = msc_species.join(meta_species, how="inner")

# Filter out species with very low abundance
df = df[(df["MscRNA"] >= 0.001) & (df["Metagenomic"] >= 0.001)]
print(f"Number of species used for correlation: {df.shape[0]}")
if df.shape[0] < 3:
    raise ValueError("Fewer than three shared species; correlation cannot be calculated.")

# Calculate Pearson correlation between the two abundance profiles
r, p = pearsonr(df["MscRNA"], df["Metagenomic"])

# Visualize the relationship between mscRNA and metagenomic abundances
plt.figure(figsize=figsize)
ax = sns.regplot(
    x=df["MscRNA"],
    y=df["Metagenomic"],
    scatter_kws=dict(
        s=point_size,
        alpha=alpha,
        edgecolor="black",
        linewidths=0.3
    ),
    line_kws=dict(color="black", linewidth=1.2),
    ci=95
)

# Label axes
ax.set_xlabel("MscRNA Abundance", fontsize=10)
ax.set_ylabel("Metagenomic Abundance", fontsize=10)
sns.despine(top=True, right=True)

# Annotate plot with correlation statistics
text = f"R={r:.2f}\nP={p:.2e}"
ax.text(
    0.05, 0.95,
    text,
    transform=ax.transAxes,
    ha="left",
    va="top",
    fontsize=10,
    bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="black", lw=0.6)
)

# Save figure as PDF
plt.tight_layout()
plt.savefig(
    "Figure3a.pdf",
    format="pdf",
    bbox_inches="tight"
)