# -*- coding: utf-8 -*-
"""
Created on Thu Aug 21 14:19:46 2025
@author: ZHENG XINGHAI
"""

import pandas as pd
import plotly.graph_objects as go
from matplotlib.colors import to_rgba

input_file, pdf_file = "cell_metadata_species.tsv", "Figure6b.pdf"
df = pd.read_csv(input_file, sep="\t")

clusters = df["cluster"].astype(str).unique().tolist()
crops = df["crop"].unique().tolist()
nodes = clusters + crops
node_indices = {x: i for i, x in enumerate(nodes)}

link_df = df.groupby(["cluster", "crop"]).size().reset_index(name="count")

cluster_color_map = {
    "0": "#A6CEE3", "1": "#98D277", "2": "#F16667",
    "3": "#FE982C", "4": "#7D54A5", "5": "#B15928"
}

crop_color_dict = {"Soybean": "#BBDED6", "Rice": "#D8BFD8", "Wheat": "#FFDAB9"}

node_colors = [cluster_color_map.get(x, "#CCCCCC") for x in clusters] + [crop_color_dict.get(x, "lightgray") for x in crops]

link_colors = [
    f"rgba({int(r*255)},{int(g*255)},{int(b*255)},0.6)"
    for r, g, b, _ in [to_rgba(cluster_color_map.get(str(x), "#CCCCCC")) for x in link_df["cluster"]]
]

fig = go.Figure(go.Sankey(
    arrangement="snap",
    node=dict(pad=50, thickness=20, line=dict(color="gray", width=0.01), label=nodes, color=node_colors),
    link=dict(
        source=link_df["cluster"].astype(str).map(node_indices),
        target=link_df["crop"].map(node_indices),
        value=link_df["count"],
        color=link_colors
    )
))

fig.update_layout(
    title=dict(text="Cluster → Crop", font=dict(size=80, family="Arial"), x=0.5),
    font=dict(size=100, family="Arial"),
    margin=dict(t=150)
)

try:
    fig.write_image(pdf_file, width=800, height=2400)
    print(f"PDF saved to: {pdf_file}")
except Exception as e:
    print("Error saving PDF. Make sure 'kaleido' is installed.")
    print(e)

fig.show()