# -*- coding: utf-8 -*-
"""
Created on Sun Sep  7 20:09:04 2025
@author: ZHENG XINGHAI
"""

import pandas as pd
import plotly.graph_objects as go
from matplotlib.colors import to_rgba

# Load cell metadata
input_file = "cell_metadata.tsv"
df = pd.read_csv(input_file, sep="\t")

# Retain the top 15 most abundant genera and group the remaining taxa as "Others"
top15_genus = df['genus'].value_counts().nlargest(15).index.tolist()
df['genus_filtered'] = df['genus'].where(df['genus'].isin(top15_genus), other="Others")

# Construct Sankey nodes representing crops, clusters, and genera
crops = df['crop'].unique().tolist()
clusters = df['cluster'].astype(str).unique().tolist()
genus = df['genus_filtered'].unique().tolist()

nodes = crops + clusters + genus
node_indices = {name: i for i, name in enumerate(nodes)}

# Calculate crop-to-cluster connections based on cell counts
link_df1 = df.groupby(['crop', 'cluster']).size().reset_index(name='count')
source1 = link_df1['crop'].map(node_indices)
target1 = link_df1['cluster'].astype(str).map(node_indices)
value1 = link_df1['count']

# Calculate cluster-to-genus connections based on cell counts
link_df2 = df.groupby(['cluster', 'genus_filtered']).size().reset_index(name='count')
source2 = link_df2['cluster'].astype(str).map(node_indices)
target2 = link_df2['genus_filtered'].map(node_indices)
value2 = link_df2['count']

# Combine all links to generate the full Sankey structure
source_indices = list(source1) + list(source2)
target_indices = list(target1) + list(target2)
values = list(value1) + list(value2)

# Define cluster color palette used to distinguish microbial transcriptional clusters
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

cluster_color_map = {c: cluster_colors.get(c, "#888888") for c in clusters}

# Define color mapping for crop types
crop_color_dict = {
    "Soybean": "#BBDED6",
    "Rice": "#D8BFD8",
    "Wheat": "#FFDAB9"
}

# Define color mapping for genera to visually distinguish microbial taxa
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

# Assign colors to nodes based on crop, cluster, and genus identity
crop_colors = [crop_color_dict.get(c, "#DDDDDD") for c in crops]
cluster_colors_list = [cluster_color_map[c] for c in clusters]
genus_colors = [genus_grouped_colors.get(g, "#E5E5E5") for g in genus]

node_colors = crop_colors + cluster_colors_list + genus_colors

# Generate link colors based on cluster identity to maintain visual consistency
link_colors = []
for s, t in zip(source_indices, target_indices):
    src_node = nodes[s]
    tgt_node = nodes[t]

    if tgt_node in clusters:
        link_colors.append(cluster_color_map[tgt_node])
    elif src_node in clusters:
        link_colors.append(cluster_color_map[src_node])
    else:
        link_colors.append("#CCCCCC")

# Apply transparency to link colors for improved visualization
link_colors = [
    f'rgba({int(r*255)},{int(g*255)},{int(b*255)},0.6)'
    for r, g, b, _ in [to_rgba(c) for c in link_colors]
]

# Italicize genus labels following standard microbiological nomenclature
genus_list = df["genus"].unique().tolist()
italic_nodes = [f"<i>{node}</i>" if node in genus_list else node for node in nodes]

# Build the Sankey diagram showing crop–cluster–genus relationships
fig = go.Figure(data=[go.Sankey(
    arrangement="snap",
    node=dict(
        pad=40,
        thickness=20,
        line=dict(color="gray", width=0.01),
        label=italic_nodes,
        color=node_colors
    ),
    link=dict(
        source=source_indices,
        target=target_indices,
        value=values,
        color=link_colors
    )
)])

# Configure figure layout and typography for publication-quality output
fig.update_layout(
    title_text="Crop → Cluster → Genus",
    title_font_size=75,
    title_font_family="Arial",
    title_x=0.5,
    font_size=50,
    font_family="Arial",
    margin=dict(t=150)
)

# Export the Sankey diagram as a high-resolution PDF figure
fig.write_image("Figure4d.pdf", width=1300, height=2000)