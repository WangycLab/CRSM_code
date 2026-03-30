# -*- coding: utf-8 -*-
"""
Created on Thu Aug 21 14:19:46 2025
@author: ZHENG XINGHAI
"""

import pandas as pd
import plotly.graph_objects as go
from matplotlib.colors import to_rgba

# Define input and output file paths
input_file = "cell_metadata_species.tsv"
pdf_file = "Figure6b.pdf"

# Load metadata table
df = pd.read_csv(input_file, sep="\t")

# Extract unique cluster and crop nodes
clusters = df['cluster'].astype(str).unique().tolist()
crops = df['crop'].unique().tolist()
nodes = clusters + crops

# Create node index mapping
node_indices = {name: i for i, name in enumerate(nodes)}

# Calculate flow counts between cluster and crop
link_df = df.groupby(['cluster', 'crop']).size().reset_index(name='count')
source_indices = link_df['cluster'].astype(str).map(node_indices)
target_indices = link_df['crop'].map(node_indices)
values = link_df['count']

# Define custom colors for each cluster
cluster_color_map = {
    "0": "#A6CEE3",
    "1": "#99CD91",
    "2": "#B89B74",
    "3": "#F06C45",
    "4": "#ED8F47",
    "5": "#825D99",
    "6": "#B15928"
}

# Define custom colors for crops
crop_color_dict = {
    "Soybean": "#BBDED6",
    "Rice": "#D8BFD8",
    "Wheat": "#FFDAB9"
}

# Assign colors to nodes: clusters with cluster colors, crops with crop colors
node_colors = [cluster_color_map.get(c, "#CCCCCC") for c in clusters] + [crop_color_dict.get(c, "lightgray") for c in crops]

# Assign colors to links based on cluster, with transparency added
link_colors = [
    f'rgba({int(r*255)},{int(g*255)},{int(b*255)},0.6)'
    for r, g, b, _ in [to_rgba(cluster_color_map.get(str(c), "#CCCCCC")) for c in link_df['cluster']]
]

# Create the Sankey diagram
fig = go.Figure(data=[go.Sankey(
    arrangement="snap",
    node=dict(
        pad=50,
        thickness=20,
        line=dict(color="gray", width=0.01),
        label=nodes,  # Add node labels here
        color=node_colors
    ),
    link=dict(
        source=source_indices,
        target=target_indices,
        value=values,
        color=link_colors
    )
)])

# Adjust layout and title appearance
fig.update_layout(
    title_text="Cluster → Crop",
    title_font_size=80,
    title_x=0.5,
    font_size=100,
    font_family="Arial",
    title_font_family="Arial",
    margin=dict(t=150)
)

# Save the figure as a PDF file
try:
    fig.write_image(pdf_file, width=800, height=2400)
    print(f"PDF saved to: {pdf_file}")
except Exception as e:
    print("Error saving PDF. Make sure 'kaleido' is installed.")
    print(e)

# Display the figure interactively
fig.show()