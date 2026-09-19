# -*- coding: utf-8 -*-
"""
Created on Thu Sep 4 16:58:41 2025
@author: ZHENG XINGHAI
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
from matplotlib.lines import Line2D

plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none", "font.family": "Arial"})

input_file, out_png, legend_png = "cell_metadata.tsv", "Figure4b_circular.pdf", "Figure4b_legend.pdf"
N_CLUSTERS, TOP_GENUS, TOP_SPECIES = 26, 15, 20
R_OUTER, R1, R2, R3, R4, R_INNER = 1.00, 0.97, 0.94, 0.91, 0.88, 0.00

cluster_annotation = {
    "0":"DNA replication|C0", "1":"DNA replication|C1", "2":"Carbohydrate transport|C2",
    "3":"Cell envelope biogenesis|C3", "4":"Carbohydrate metabolism|C4", "5":"Phosphorylation|C5",
    "6":"DNA replication|C6", "7":"Cell envelope biogenesis|C7", "8":"Signal transduction|C8",
    "9":"Antibiotic resistance|C9", "10":"Protein turnover|C10", "11":"Stress response|C11",
    "12":"Protein insertion|C12", "13":"Cell wall remodeling|C13", "14":"Respiration|C14",
    "15":"Protein folding|C15", "16":"DNA repair|C16", "17":"DNA replication|C17",
    "18":"DNA segregation|C18", "19":"Carbohydrate metabolism|C19", "20":"Protein folding|C20",
    "21":"DNA segregation|C21", "22":"Stress response|C22", "23":"DNA segregation|C23",
    "24":"Protein degradation|C24", "25":"Phosphate transport|C25"
}

cluster_colors = {
    "0":"#A6CEE3", "1":"#6AA8CE", "2":"#2F82B9", "3":"#4E98A6", "4":"#8EC694",
    "5":"#98D277", "6":"#60B64D", "7":"#439F34", "8":"#9B9C64", "9":"#F29A94",
    "10":"#F16667", "11":"#E62E30", "12":"#EA4833", "13":"#F59057", "14":"#FDB45D",
    "15":"#FE982C", "16":"#FC8108", "17":"#E59766", "18":"#CEADC4", "19":"#A787C0",
    "20":"#7D54A5", "21":"#8D6B99", "22":"#CFC099", "23":"#F5EB8B", "24":"#D3A259",
    "25":"#B15928"
}

crop_colors = {"Soybean":"#BBDED6", "Rice":"#D8BFD8", "Wheat":"#FFDAB9"}

genus_grouped_colors = {
    "Sinorhizobium":"#1B9E77", "Pseudomonas":"#7A7E3C", "Burkholderia":"#D95F02",
    "Brucella":"#A7675A", "Pedobacter":"#7570B3", "Salmonella":"#AD4C9E",
    "Escherichia":"#E7298A", "Enhydrobacter":"#A66753", "Rhizobium":"#66A61E",
    "Duganella":"#A6A810", "Massilia":"#E6AB02", "Agrobacterium":"#C5900F",
    "Methylococcus":"#A6761D", "Bacillus":"#866E41", "Nitrosomonas":"#666666",
    "Others":"#E5E5E5"
}

species_grouped_colors = {
    "Sinorhizobium fredii":"#8DD3C7", "Salmonella enterica":"#CFECBB",
    "Pseudomonas aeruginosa":"#F4F4B9", "Escherichia coli":"#CFCCCF",
    "Brucella melitensis":"#D1A7B9", "Burkholderia pseudomallei":"#F4867C",
    "Pedobacter sp. FW305-3-2-15-E-R2A2":"#C0979F", "Enhydrobacter sp.":"#86B1CD",
    "Pseudomonas chlororaphis":"#CEB28B", "Burkholderia cepacia":"#EDBC63",
    "Burkholderia thailandensis":"#C2D567", "Brucella abortus":"#CDD796",
    "Aquirufa nivalisilvae":"#F8CDDE", "Methylococcus sp. EFPC2":"#E9D3DE",
    "Duganella zoogloeoides":"#D5CFD6", "Sinorhizobium meliloti":"#C59CC5",
    "Rugamonas sp. DEMB1":"#C09CBF", "Agrobacterium tumefaciens":"#C9DAC3",
    "Emticicia sp. 21SJ11W-3":"#E1EBA0", "Pseudomonas sp. S150":"#FFED6F",
    "Others":"#E5E5E5"
}

df = pd.read_csv(input_file, sep="\t")
df["cluster"] = df["cluster"].astype(str)
df = df[df["cluster"].isin(cluster_annotation)].copy()
df["cluster_annotation"] = df["cluster"].map(cluster_annotation)

top_clusters = [c for c in cluster_annotation if c in df["cluster"].unique()]
annotation_order = [cluster_annotation[c] for c in top_clusters]
crops = df["crop"].astype(str).unique().tolist()

top_genus_list = df["genus"].value_counts().nlargest(TOP_GENUS).index.tolist()
top_species_list = df["species"].value_counts().nlargest(TOP_SPECIES).index.tolist()

df["genus_filtered"] = np.where(df["genus"].isin(top_genus_list), df["genus"], "Others")
df["species_filtered"] = np.where(df["species"].isin(top_species_list), df["species"], "Others")

cluster_to_crop, cluster_to_genus, cluster_to_species = {}, {}, {}

for c in top_clusters:
    sub = df[df["cluster"] == c]
    cluster_to_crop[c] = sub["crop"].value_counts(normalize=True).reindex(crops, fill_value=0).to_dict()
    cluster_to_genus[c] = sub["genus_filtered"].value_counts(normalize=True).reindex(top_genus_list + ["Others"], fill_value=0).to_dict()
    cluster_to_species[c] = sub["species_filtered"].value_counts(normalize=True).reindex(top_species_list + ["Others"], fill_value=0).to_dict()

def draw_ring(ax, start_deg, end_deg, r_inner, r_outer, parts_dict, color_map):
    current = start_deg
    for key, frac in parts_dict.items():
        if frac <= 0: continue
        theta2 = current + frac * (end_deg - start_deg)
        ax.add_patch(Wedge((0, 0), r_outer, current, theta2, width=r_outer-r_inner,
                           facecolor=color_map.get(key, "grey"), edgecolor="white", linewidth=0.5))
        current = theta2

fig, ax = plt.subplots(figsize=(12, 12), subplot_kw={"aspect": "equal"})
ax.set_xlim(-1.05, 1.05); ax.set_ylim(-1.05, 1.05); ax.axis("off")

angle_per_cluster = 360 / len(top_clusters)

for i, c in enumerate(top_clusters):
    theta1, theta2 = i * angle_per_cluster, (i + 1) * angle_per_cluster
    ax.add_patch(Wedge((0, 0), R_OUTER, theta1, theta2, width=R_OUTER-R1,
                       facecolor=cluster_colors.get(c, "grey"), edgecolor="white", linewidth=0.8))
    draw_ring(ax, theta1, theta2, R1, R2, cluster_to_crop[c], crop_colors)
    draw_ring(ax, theta1, theta2, R2, R3, cluster_to_genus[c], genus_grouped_colors)
    draw_ring(ax, theta1, theta2, R3, R4, cluster_to_species[c], species_grouped_colors)

if R_INNER > 0:
    ax.add_patch(Wedge((0, 0), R4, 0, 360, width=R4-R_INNER, facecolor="white", edgecolor="white"))

plt.savefig(out_png, dpi=300, bbox_inches="tight")
plt.close(fig)

def make_italic_label(s): return s

legend_fig, legend_ax = plt.subplots(figsize=(15, 12))
legend_ax.axis("off")

legends = [
    ("Functional Cluster", [
        Line2D([0], [0], marker="o", linestyle="", markersize=8,
               markerfacecolor=cluster_colors.get(c, "grey"), markeredgecolor="none",
               label=cluster_annotation[c]) for c in top_clusters
    ], 4),
    ("Crop", [
        Line2D([0], [0], marker="o", linestyle="", markersize=8,
               markerfacecolor=crop_colors.get(c, "grey"), markeredgecolor="none", label=c)
        for c in crops
    ], 3),
    ("Genus", [
        Line2D([0], [0], marker="o", linestyle="", markersize=8,
               markerfacecolor=genus_grouped_colors.get(g, "grey"), markeredgecolor="none",
               label=make_italic_label(g)) for g in top_genus_list + ["Others"]
    ], 4),
    ("Species", [
        Line2D([0], [0], marker="o", linestyle="", markersize=8,
               markerfacecolor=species_grouped_colors.get(s, "grey"), markeredgecolor="none",
               label=make_italic_label(s)) for s in top_species_list + ["Others"]
    ], 3)
]

y0 = 1.0
for title, handles, ncol in legends:
    leg = legend_ax.legend(handles=handles, title=title, loc="upper center",
                           bbox_to_anchor=(0.5, y0), ncol=ncol, frameon=False,
                           fontsize=12, title_fontsize=20)
    if title in {"Genus", "Species"}:
        for text in leg.get_texts(): text.set_fontstyle("italic")
    legend_ax.add_artist(leg)
    y0 -= {"Functional Cluster": 0.25, "Crop": 0.10, "Genus": 0.15, "Species": 0.20}[title]

legend_fig.savefig(legend_png, dpi=300, bbox_inches="tight")
plt.close(legend_fig)