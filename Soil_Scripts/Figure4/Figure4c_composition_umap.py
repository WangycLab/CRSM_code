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

# Configure matplotlib to keep text editable in vector graphics software
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams["font.family"] = "Arial"

# Define input metadata file and output figure path
input_file = "cell_metadata.tsv"
out_png = "Figure4c_circular.pdf"

# Set parameters controlling the number of displayed taxa and ring radii
N_CLUSTERS = 20
TOP_GENUS = 15
TOP_SPECIES = 20
R_OUTER = 1.00
R1 = 0.97
R2 = 0.94
R3 = 0.91
R4 = 0.88
R_INNER = 0.00

plt.rcParams["font.family"] = "Arial"

# Define color palette for transcriptional clusters
cluster_colors = {
    "0": "#A6CEE3",
    "1": "#579CC7",
    "2": "#3688AD",
    "3": "#8BC395",
    "4": "#89CB6C",
    "5": "#40A635",
    "6": "#919D5F",
    "7": "#F99392",
    "8": "#EB494A",
    "9": "#E83C2D",
    "10": "#F79C5D",
    "11": "#FDA746",
    "12": "#FE8205",
    "13": "#E39970",
    "14": "#BFA5CF",
    "15": "#8861AC",
    "16": "#917099",
    "17": "#E7E099",
    "18": "#DEB969",
    "19": "#B15928"
}

# Define color palette for crop categories
crop_colors = {
    "Soybean": "#BBDED6",
    "Rice": "#D8BFD8",
    "Wheat": "#FFDAB9"
}

# Define color palette for dominant microbial genera
genus_grouped_colors = {
    "Sinorhizobium": "#1B9E77",
    "Pseudomonas": "#7A7E3C",
    "Burkholderia": "#D95F02",
    "Brucella": "#A7675A",
    "Pedobacter": "#7570B3",
    "Salmonella": "#AD4C9E",
    "Escherichia": "#E7298A",
    "Enhydrobacter": "#A66753",
    "Rhizobium": "#66A61E",
    "Massilia": "#A6A810",
    "Duganella": "#E6AB02",
    "Agrobacterium": "#C5900F",
    "Nitrosomonas": "#A6761D",
    "Methylococcus": "#866E41",
    "Bacillus": "#666666",
    "Others": "#E5E5E5"
}

# Define color palette for dominant microbial species
species_grouped_colors = {
    "Sinorhizobium fredii": "#8DD3C7",
    "Salmonella enterica": "#CFECBB",
    "Pseudomonas aeruginosa": "#F4F4B9",
    "Escherichia coli": "#CFCCCF",
    "Brucella melitensis": "#D1A7B9",
    "Burkholderia pseudomallei": "#F4867C",
    "Pedobacter sp. FW305": "#C0979F",
    "Enhydrobacter sp.": "#86B1CD",
    "Pseudomonas chlororaphis": "#CEB28B",
    "Burkholderia cepacia": "#EDBC63",
    "Burkholderia thailandensis": "#C2D567",
    "Brucella abortus": "#CDD796",
    "Aquirufa nivalisilvae": "#F8CDDE",
    "Duganella zoogloeoides": "#E9D3DE",
    "Methylococcus sp. EFPC2": "#D5CFD6",
    "Sinorhizobium meliloti": "#C59CC5",
    "Agrobacterium tumefaciens": "#C09CBF",
    "Rugamonas sp. DEMB1": "#C9DAC3",
    "Pseudomonas sp. S150": "#E1EBA0",
    "Emticicia sp. 21SJ11W-3": "#FFED6F",
    "Others": "#E5E5E5"
}

# Load cell metadata containing cluster, crop, genus, and species annotations
df = pd.read_csv(input_file, sep="\t")

# Select the most abundant clusters and convert cluster identifiers to strings
cluster_counts = df["cluster"].value_counts()
top_clusters = cluster_counts.nlargest(N_CLUSTERS).index.astype(str).tolist()
df = df[df["cluster"].astype(str).isin(top_clusters)].copy()
df["cluster"] = df["cluster"].astype(str)

# Extract crop categories present in the dataset
crops = df["crop"].astype(str).unique().tolist()

# Retain the most abundant genera and species while grouping others into an "Others" category
top_genus_list = df["genus"].value_counts().nlargest(TOP_GENUS).index.tolist()
df["genus_filtered"] = np.where(df["genus"].isin(top_genus_list), df["genus"], "Others")
top_species_list = df["species"].value_counts().nlargest(TOP_SPECIES).index.tolist()
df["species_filtered"] = np.where(df["species"].isin(top_species_list), df["species"], "Others")

# Calculate proportional composition of crop, genus, and species within each cluster
cluster_to_crop = {}
cluster_to_genus = {}
cluster_to_species = {}
for c in top_clusters:
    sub = df[df["cluster"] == c]
    cluster_to_crop[c] = sub["crop"].value_counts(normalize=True).reindex(crops, fill_value=0).to_dict()
    cluster_to_genus[c] = sub["genus_filtered"].value_counts(normalize=True).reindex(list(top_genus_list)+["Others"], fill_value=0).to_dict()
    cluster_to_species[c] = sub["species_filtered"].value_counts(normalize=True).reindex(list(top_species_list)+["Others"], fill_value=0).to_dict()

# Define a helper function to draw segmented circular ring layers
def draw_ring(ax, start_deg, end_deg, r_inner, r_outer, parts_dict, color_map):
    angle_span = end_deg - start_deg
    current = start_deg
    for k, frac in parts_dict.items():
        if frac <= 0:
            continue
        theta1 = current
        theta2 = current + frac * angle_span
        wedge = Wedge(center=(0, 0), r=r_outer, theta1=theta1, theta2=theta2,
                      width=(r_outer - r_inner), facecolor=color_map.get(k, "grey"),
                      edgecolor="white", linewidth=0.5)
        ax.add_patch(wedge)
        current = theta2

# Create circular plot canvas for multi-layer cluster composition visualization
fig, ax = plt.subplots(figsize=(12, 12), subplot_kw=dict(aspect="equal"))
ax.set_xlim(-1.05, 1.05)
ax.set_ylim(-1.05, 1.05)
ax.axis("off")

# Determine angular width assigned to each cluster segment
n = len(top_clusters)
angle_per_cluster = 360.0 / n

# Draw the outermost ring representing cluster identities
for i, c in enumerate(top_clusters):
    theta1 = i * angle_per_cluster
    theta2 = (i + 1) * angle_per_cluster
    wedge = Wedge(center=(0, 0), r=R_OUTER, theta1=theta1, theta2=theta2,
                  width=(R_OUTER - R1), facecolor=cluster_colors.get(c, "grey"),
                  edgecolor="white", linewidth=0.8)
    ax.add_patch(wedge)

# Draw inner rings representing crop composition, genus composition, and species composition
for i, c in enumerate(top_clusters):
    theta1 = i * angle_per_cluster
    theta2 = (i + 1) * angle_per_cluster
    draw_ring(ax, theta1, theta2, R1, R2, cluster_to_crop[c], crop_colors)
    draw_ring(ax, theta1, theta2, R2, R3, cluster_to_genus[c], genus_grouped_colors)
    draw_ring(ax, theta1, theta2, R3, R4, cluster_to_species[c], species_grouped_colors)

# Optionally draw a central core region if an inner radius is specified
if R_INNER > 0:
    core = Wedge(center=(0, 0), r=R4, theta1=0, theta2=360,
                 width=R4 - R_INNER, facecolor="white", edgecolor="white")
    ax.add_patch(core)

# Export the circular composition plot as a publication-quality PDF
plt.savefig(out_png, dpi=300, bbox_inches="tight")
plt.close(fig)

# Define helper function for formatting italic labels in legends
def make_italic_label(s):
    return s

# Create a separate figure containing legends for clusters, crops, genera, and species
legend_fig, legend_ax = plt.subplots(figsize=(15, 12))
legend_ax.axis("off")

# Prepare legend handles for clusters, crops, genera, and species
legends = []
cluster_handles = [Line2D([0], [0], marker='o', linestyle='', markersize=8,
                          markerfacecolor=cluster_colors.get(c, "grey"),
                          markeredgecolor='none', label=f"Cluster {c}") for c in top_clusters]
legends.append(("Cluster", cluster_handles, 5))

crop_handles = [Line2D([0], [0], marker='o', linestyle='', markersize=8,
                       markerfacecolor=crop_colors.get(h, "grey"),
                       markeredgecolor='none', label=h) for h in crops]
legends.append(("Crop", crop_handles, 3))

genus_handles = [Line2D([0], [0], marker='o', linestyle='', markersize=8,
                        markerfacecolor=genus_grouped_colors.get(g, "grey"),
                        markeredgecolor='none', label=make_italic_label(g)) for g in list(top_genus_list)+["Others"]]
legends.append(("Genus", genus_handles, 4))

species_handles = [Line2D([0], [0], marker='o', linestyle='', markersize=8,
                          markerfacecolor=species_grouped_colors.get(s, "grey"),
                          markeredgecolor='none', label=make_italic_label(s)) for s in list(top_species_list)+["Others"]]
legends.append(("Species", species_handles, 3))

# Arrange legends vertically and apply italic formatting to genus and species labels
y0 = 1.0
for title, handles, ncol in legends:
    leg = legend_ax.legend(handles=handles, title=title, loc="upper center",
                           bbox_to_anchor=(0.5, y0), ncol=ncol, frameon=False,
                           fontsize=12, title_fontsize=20)

    if title in ["Genus", "Species"]:
        for text in leg.get_texts():
            text.set_fontstyle("italic")

    legend_ax.add_artist(leg)

    if title == "Crop":
        y0 -= 0.1
    elif title == "Cluster":
        y0 -= 0.15
    elif title == "Genus":
        y0 -= 0.16
    else:
        y0 -= 0.2

# Save the legend panel as a separate high-resolution PDF
legend_png = "Figure4c_legend.pdf"
legend_fig.savefig(legend_png, dpi=300, bbox_inches="tight")
plt.close(legend_fig)