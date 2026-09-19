# -*- coding: utf-8 -*-
"""
Created on Sun Sep  7 20:09:04 2025
@author: ZHENG XINGHAI
"""

import pandas as pd
import plotly.graph_objects as go
from matplotlib.colors import to_rgba

df = pd.read_csv("cell_metadata.tsv", sep="\t")

top15_genus = df["genus"].value_counts().nlargest(15).index.tolist()
df["genus_filtered"] = df["genus"].where(df["genus"].isin(top15_genus), "Others")

crops = df["crop"].unique().tolist()
clusters = df["cluster"].astype(str).unique().tolist()
genus = df["genus_filtered"].unique().tolist()
nodes = crops + clusters + genus
node_indices = {n: i for i, n in enumerate(nodes)}

link_df1 = df.groupby(["crop", "cluster"]).size().reset_index(name="count")
link_df2 = df.groupby(["cluster", "genus_filtered"]).size().reset_index(name="count")

source = list(link_df1["crop"].map(node_indices)) + list(link_df2["cluster"].astype(str).map(node_indices))
target = list(link_df1["cluster"].astype(str).map(node_indices)) + list(link_df2["genus_filtered"].map(node_indices))
values = list(link_df1["count"]) + list(link_df2["count"])

cluster_colors = {
    "0":"#A6CEE3", "1":"#6AA8CE", "2":"#2F82B9", "3":"#4E98A6", "4":"#8EC694",
    "5":"#98D277", "6":"#60B64D", "7":"#439F34", "8":"#9B9C64", "9":"#F29A94",
    "10":"#F16667", "11":"#E62E30", "12":"#EA4833", "13":"#F59057", "14":"#FDB45D",
    "15":"#FE982C", "16":"#FC8108", "17":"#E59766", "18":"#CEADC4", "19":"#A787C0",
    "20":"#7D54A5", "21":"#8D6B99", "22":"#CFC099", "23":"#F5EB8B", "24":"#D3A259",
    "25":"#B15928"
}

crop_colors = {"Soybean":"#BBDED6", "Rice":"#D8BFD8", "Wheat":"#FFDAB9"}

genus_colors = {
    "Sinorhizobium":"#1B9E77", "Pseudomonas":"#7A7E3C", "Burkholderia":"#D95F02",
    "Brucella":"#A7675A", "Pedobacter":"#7570B3", "Salmonella":"#AD4C9E",
    "Escherichia":"#E7298A", "Enhydrobacter":"#A66753", "Rhizobium":"#66A61E",
    "Duganella":"#A6A810", "Massilia":"#E6AB02", "Agrobacterium":"#C5900F",
    "Methylococcus":"#A6761D", "Bacillus":"#866E41", "Nitrosomonas":"#666666",
    "Others":"#E5E5E5"
}

cluster_color_map = {c: cluster_colors.get(c, "#888888") for c in clusters}
node_colors = (
    [crop_colors.get(c, "#DDDDDD") for c in crops] +
    [cluster_color_map[c] for c in clusters] +
    [genus_colors.get(g, "#E5E5E5") for g in genus]
)

link_colors = []
for s, t in zip(source, target):
    src, tgt = nodes[s], nodes[t]
    link_colors.append(cluster_color_map.get(tgt, cluster_color_map.get(src, "#CCCCCC")))

link_colors = [
    f"rgba({int(r*255)},{int(g*255)},{int(b*255)},0.6)"
    for r, g, b, _ in map(to_rgba, link_colors)
]

genus_list = set(df["genus"].astype(str))
node_labels = [f"<i>{n}</i>" if n in genus_list else n for n in nodes]

fig = go.Figure(go.Sankey(
    arrangement="snap",
    node=dict(
        pad=40, thickness=20, line=dict(color="gray", width=0.01),
        label=node_labels, color=node_colors
    ),
    link=dict(source=source, target=target, value=values, color=link_colors)
))

fig.update_layout(
    title_text="Crop → Cluster → Genus", title_font_size=75, title_font_family="Arial", title_x=0.5,
    font_size=50, font_family="Arial", margin=dict(t=150)
)

fig.write_image("Figure4c.pdf", width=1300, height=2000)