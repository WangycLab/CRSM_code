# -*- coding: utf-8 -*-
"""
Created on Mon Sep 29 18:53:18 2025
@author: ZHENG XINGHAI
"""

import pandas as pd
import os, glob
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none", "font.family": "Arial"})

AXIS_LABEL_SIZE, TICK_LABEL_SIZE = 25, 18
CBAR_WIDTH, CBAR_LENGTH, CBAR_TICK_SIZE, CBAR_LABEL_SIZE = 0.1, 1, 15, 18
module_order = ["blue", "brown", "yellow", "red", "black", "turquoise", "green"]

input_dir = "COG_enrichment_by_module"
files = glob.glob(os.path.join(input_dir, "COG_enrichment_module_*.csv"))

p_value_dict, rich_factor_dict = {}, {}
for file in files:
    module = os.path.basename(file).replace("COG_enrichment_module_", "").replace(".csv", "")
    df = pd.read_csv(file)
    df = df[(df["overlap"] > 0) & (df["rich_factor"] > 0.1)]
    p_value_dict[module] = df.set_index("COG_ID")["p_value"]
    rich_factor_dict[module] = df.set_index("COG_ID")["rich_factor"]

p_value_df = pd.DataFrame(p_value_dict).fillna(1).reindex(columns=module_order)
rich_factor_df = pd.DataFrame(rich_factor_dict).fillna(0).reindex(columns=module_order)

cog_anno = pd.read_csv("cog-24.def.tab", sep="\t", header=None, usecols=[0, 2], on_bad_lines="skip", engine="python")
cog_anno.columns = ["COG_ID", "Annotation"]
cog_anno_dict = cog_anno.set_index("COG_ID")["Annotation"].to_dict()
rename_cog = lambda x: cog_anno_dict.get(x, x).strip()
rich_factor_df.index, p_value_df.index = rich_factor_df.index.map(rename_cog), p_value_df.index.map(rename_cog)

rich_factor_df_norm = rich_factor_df.copy()
for col in rich_factor_df_norm:
    mn, mx = rich_factor_df_norm[col].min(), rich_factor_df_norm[col].max()
    rich_factor_df_norm[col] = (rich_factor_df_norm[col] - mn) / (mx - mn) if mx > mn else 0

cg = sns.clustermap(rich_factor_df_norm, method="average", metric="euclidean", row_cluster=True, col_cluster=False, cmap="Blues", figsize=(1, 1))
row_order = cg.dendrogram_row.reordered_ind
plt.close()
rich_factor_df_norm, p_value_df = rich_factor_df_norm.iloc[row_order], p_value_df.iloc[row_order]

abbreviation_dict = {
    'Cytoplasmic potassium-binding protein Kbp/XkdP/YgaU, contains LysM domain': 'Kbp potassium protein',
    'Phage tail tape-measure protein, controls tail length': 'Phage tail protein',
    'Phage portal protein BeeE': 'Phage portal protein',
    'HrpA-like RNA helicase': 'HrpA RNA helicase',
    'Transposase, IS1182 family': 'IS1182 transposase',
    'tRNA-C32 2-thiocytidine or tRNA(Ile)-C34 C2-lysylcytidine synthase TtcA/TilS/MesJ': 'tRNA modification enzyme',
    'Outer membrane receptor for ferric coprogen and ferric-rhodotorulic acid': 'Ferric coprogen receptor',
    'FMN-dependent dehydrogenase, includes L-lactate dehydrogenase and type II isopentenyl diphosphate isomerase': 'FMN dehydrogenase',
    'Rhodanese-related sulfurtransferase': 'Rhodanese sulfurtransferase',
    'Mg2+ efflux pump MpfA, contains CBS pair and CorC-HlyC domains, TlyC/UPF0053 family': 'Mg2+ efflux pump',
    'NADH dehydrogenase, FAD-containing subunit': 'NADH dehydrogenase',
    'Outer membrane usher protein FimD/PapC': 'Outer membrane usher',
    'Zinc transporter ZupT': 'Zinc transporter',
    'Glutaredoxin': 'Glutaredoxin',
    'HNH family endonuclease, includes 5-methylcytosine-specific restriction endonuclease McrA': 'HNH endonuclease',
    'Formylglycine-generating enzyme, required for sulfatase activity, contains SUMF1/FGE domain': 'Formylglycine enzyme',
    'Outer membrane protein TolC': 'Outer membrane protein',
    'Signal recognition particle GTPase FtsY': 'SRP GTPase FtsY',
    'Cbb3-type cytochrome oxidase, cytochrome c subunit FixO': 'Cytochrome oxidase FixO',
    'ParA-like ATPase involved in chromosome/plasmid partitioning or cellulose biosynthesis protein BcsQ': 'ParA-like ATPase',
    'ABC-type oligopeptide transport system, periplasmic component': 'ABC oligopeptide transporter',
    'Small heat shock protein IbpA, HSP20 family': 'Small heat shock protein',
    'Uncharacterized conserved protein': 'Conserved protein',
    'RecA/RadA recombinase': 'RecA recombinase',
    'Integral membrane protein YbhL, putative Ca2+ regulator, Bax inhibitor (BI-1)/TMBIM family': 'YbhL membrane protein',
    'Membrane protein insertase Oxa1/YidC/SpoIIIJ': 'Membrane protein insertase',
    'Predicted lipid-binding transport protein, Tim44 family': 'Tim44 transport protein',
    'Phosphoglycerate dehydrogenase or related dehydrogenase': 'Phosphoglycerate dehydrogenase',
    'Zn-dependent oligopeptidase, M3 family': 'Zn oligopeptidase',
    'ABC-type branched-chain amino acid transport system, permease component': 'ABC BCAA permease',
    'Superfamily II DNA or RNA helicase, SNF2 family': 'SNF2 helicase',
    'K+ uptake protein Kup': 'K+ uptake protein',
    'HD superfamily phosphohydrolase': 'HD phosphohydrolase',
    'Phosphopantetheinyl transferase': 'Phosphopantetheinyl transferase'
}

shorten_label = lambda x: next((v for k, v in abbreviation_dict.items() if k in x), x)
rich_factor_df_norm.index, p_value_df.index = rich_factor_df_norm.index.map(shorten_label), p_value_df.index.map(shorten_label)

plt.figure(figsize=(7.5, 12))
ax = sns.heatmap(rich_factor_df_norm, cmap="Purples", linewidths=1, linecolor="grey", annot=False,
                 cbar_kws={"label": "Rich Factor", "fraction": CBAR_WIDTH, "shrink": CBAR_LENGTH})

for i, j in zip(*((p_value_df < 0.05).to_numpy()).nonzero()):
    ax.text(j + 0.5, i + 0.95, "*", ha="center", va="center", color="white", fontsize=30, fontweight="bold")

ax.set(xlabel="Module", ylabel="COG Function")
ax.tick_params(axis="x", labelsize=TICK_LABEL_SIZE, rotation=90)
ax.tick_params(axis="y", labelsize=TICK_LABEL_SIZE)

cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=CBAR_TICK_SIZE)
cbar.set_label("Scaled Rich Factor", fontsize=CBAR_LABEL_SIZE)

plt.tight_layout()
plt.savefig("Figure5e.pdf", dpi=300, bbox_inches="tight")
plt.show()