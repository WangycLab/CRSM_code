# -*- coding: utf-8 -*-
"""
Created on Fri Jul 25 16:22:45 2025
@author: ZHENG XINGHAI & XIONG XIAO
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

# Settings
plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.size"] = 9
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
sns.set(style="ticks", context="paper")
figsize = (3.5, 3.8)
cutoffs = [0.001]

# Load data
msc = pd.read_csv("msc_species.tsv", sep="\t")
meta = pd.read_csv("meta_speices.tsv", sep="\t")

def calculate_metrics(cutoff):
    msc_df = msc.loc[msc["fraction_total_reads"] >= cutoff, ["name", "fraction_total_reads"]]
    meta_df = meta.loc[meta["fraction_total_reads"] >= cutoff, ["name", "fraction_total_reads"]]
    msc_set, meta_set = set(msc_df["name"]), set(meta_df["name"])
    shared, union = msc_set & meta_set, msc_set | meta_set

    jaccard = len(shared) / len(union) if union else np.nan
    recall_meta = len(shared) / len(meta_set) if meta_set else np.nan
    recall_msc = len(shared) / len(msc_set) if msc_set else np.nan

    msc_ab = msc_df.set_index("name").rename(columns={"fraction_total_reads": "MscRNA"})
    meta_ab = meta_df.set_index("name").rename(columns={"fraction_total_reads": "Metagenomic"})
    df = msc_ab.join(meta_ab, how="inner")
    df["log10_MscRNA"], df["log10_Meta"] = np.log10(df["MscRNA"]), np.log10(df["Metagenomic"])

    r, p = pearsonr(df["log10_MscRNA"], df["log10_Meta"]) if len(df) >= 3 else (np.nan, np.nan)

    return {
        "cutoff": cutoff,
        "mscRNA species": len(msc_set),
        "Meta species": len(meta_set),
        "Shared species": len(shared),
        "Jaccard": jaccard,
        "Recall Meta→mscRNA": recall_meta,
        "Recall mscRNA→Meta": recall_msc,
        "Pearson R": r,
        "P": p,
        "df": df
    }

# Statistics
results = [calculate_metrics(c) for c in cutoffs]
summary_table = pd.DataFrame(results).drop(columns="df")
print("\n============================\nCutoff statistics\n============================")
print(summary_table.to_string(index=False))
summary_table.to_csv("cutoff_statistics.tsv", sep="\t", index=False)

# Plot
res, df = results[0], results[0]["df"]
fig, ax = plt.subplots(figsize=figsize)

ax.scatter(df["log10_MscRNA"], df["log10_Meta"], s=35, alpha=0.8, edgecolors="black", linewidths=0.3)

x, y = df["log10_MscRNA"].values, df["log10_Meta"].values
if len(x) >= 2:
    slope, intercept = np.polyfit(x, y, 1)
    x_line = np.linspace(x.min(), x.max(), 200)
    ax.plot(x_line, slope * x_line + intercept, color="black", linewidth=1)

ax.set_xlabel("log10 mscRNA species relative abundance", fontsize=9)
ax.set_ylabel("log10 metagenomic species relative abundance", fontsize=9)

text = (
    f"n={res['Shared species']}\n"
    f"R={res['Pearson R']:.2f}\n"
    f"P={res['P']:.2e}\n\n"
    f"Jaccard={res['Jaccard']:.2f}\n"
    f"Recall (Meta→mscRNA)={res['Recall Meta→mscRNA']:.2f}\n"
    f"Recall (mscRNA→Meta)={res['Recall mscRNA→Meta']:.2f}"
)

ax.text(0.05, 0.95, text, transform=ax.transAxes, fontsize=7, va="top", ha="left",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="black", linewidth=0.5))

sns.despine(ax=ax)
plt.tight_layout()
plt.savefig("Figure3a.pdf", dpi=600, bbox_inches="tight")
plt.show()