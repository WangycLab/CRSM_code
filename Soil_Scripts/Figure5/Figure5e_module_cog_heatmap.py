# -*- coding: utf-8 -*-
"""
Created on Mon Sep 29 18:53:18 2025
@author: ZHENG XINGHAI
"""

import pandas as pd
import os
import glob
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

# Configure matplotlib to ensure editable
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams['font.family'] = 'Arial'

# Define plotting parameters and module order
AXIS_LABEL_SIZE = 25
TICK_LABEL_SIZE = 18
CBAR_WIDTH = 0.1
CBAR_LENGTH = 1
CBAR_TICK_SIZE = 15
CBAR_LABEL_SIZE = 18
module_order = ["blue", "brown", "green", "red", "turquoise", "yellow"]

# Read COG enrichment results for all modules
input_dir = "COG_enrichment_by_module"
files = glob.glob(os.path.join(input_dir, "COG_enrichment_module_*.csv"))

# Filter results and organize p-values and rich factors by module
p_value_dict = {}
rich_factor_dict = {}
for file in files:
    module_name = os.path.basename(file).replace("COG_enrichment_module_", "").replace(".csv", "")
    df = pd.read_csv(file)
    df = df[(df["overlap"] >= 1) & (df["rich_factor"] >= 0.3)]
    p_value_dict[module_name] = df.set_index("COG_ID")["p_value"]
    rich_factor_dict[module_name] = df.set_index("COG_ID")["rich_factor"]

# Merge enrichment data into matrices with modules as columns
p_value_df = pd.DataFrame(p_value_dict).fillna(1)
rich_factor_df = pd.DataFrame(rich_factor_dict).fillna(0)
rich_factor_df = rich_factor_df.reindex(columns=module_order)
p_value_df = p_value_df.reindex(columns=module_order)

# Load COG annotation and map COG IDs to descriptive function names
cog_anno = pd.read_csv(
    "cog-24.def.tab",
    sep="\t",
    header=None,
    usecols=[0, 2],
    on_bad_lines="skip",
    engine="python"
)
cog_anno.columns = ["COG_ID", "Annotation"]
cog_anno_dict = cog_anno.set_index("COG_ID")["Annotation"].to_dict()
rich_factor_df.index = rich_factor_df.index.map(lambda x: cog_anno_dict.get(x, x)).str.strip()
p_value_df.index = p_value_df.index.map(lambda x: cog_anno_dict.get(x, x)).str.strip()

# Normalize rich factor values column-wise to 0-1 scale
rich_factor_df_norm = rich_factor_df.copy()
for col in rich_factor_df_norm.columns:
    col_min = rich_factor_df_norm[col].min()
    col_max = rich_factor_df_norm[col].max()
    if col_max > col_min:
        rich_factor_df_norm[col] = (rich_factor_df_norm[col] - col_min) / (col_max - col_min)
    else:
        rich_factor_df_norm[col] = 0

# Cluster rows based on rich factor patterns to determine plotting order
cg = sns.clustermap(
    rich_factor_df_norm,
    method="average",
    metric="euclidean",
    row_cluster=True,
    col_cluster=False,
    cmap="Blues",
    figsize=(1, 1)
)
row_order = cg.dendrogram_row.reordered_ind
plt.close()
rich_factor_df_norm = rich_factor_df_norm.iloc[row_order, :]
p_value_df = p_value_df.iloc[row_order, :]

# Define a dictionary to shorten long COG function names for y-axis
abbreviation_dict = {
'HrpA-like RNA helicase': 'HrpA RNA helicase',
'Transposase, IS1182 family': 'Transposase IS1182',
'Cytoplasmic potassium-binding protein Kbp/XkdP/YgaU, contains LysM domain': 'Kbp potassium-binding protein',
'Periplasmic deferrochelatase/peroxidase EfeB': 'Periplasmic EfeB',
'Alcohol dehydrogenase YqhD, Fe-dependent ADH family': 'Alcohol dehydrogenase YqhD',
'Transposase': 'Transposase',
'Nitric oxide response protein NnrS': 'NO response protein NnrS',
'Phosphoribosylcarboxyaminoimidazole (NCAIR) mutase': 'NCAIR mutase',
'Inorganic triphosphatase YgiF, contains CYTH and CHAD domains': 'Triphosphatase YgiF',
'FoF1-type ATP synthase, membrane subunit c/Archaeal/vacuolar-type H+-ATPase, subunit K': 'ATP synthase subunit c',
'5-oxoprolinase subunit C/Allophanate hydrolase subunit 2': 'Oxoprolinase subunit C',
'Uncharacterized secreted protein, contains LVIVD repeats, choice-of-anchor domain': 'Secreted LVIVD protein',
'tRNA-C32 2-thiocytidine or tRNA(Ile)-C34 C2-lysylcytidine synthase TtcA/TilS/MesJ': 'tRNA modification enzyme',
'DNA-binding response regulator, LytR/AlgR family': 'Response regulator LytR',
'Phage portal protein BeeE': 'Phage portal BeeE',
'ABC-type arginine/histidine transport system, permease component': 'ABC Arg/His permease',
'Phage tail tape-measure protein, controls tail length': 'Phage tail tape protein',
'Fructose-bisphosphate aldolase class 1': 'Fructose-bisphosphate aldolase',
'Gluconate kinase': 'Gluconate kinase',
'Signal recognition particle GTPase FtsY': 'SRP GTPase FtsY',
"tRNA(Leu) C34 or U34 (ribose-2'-O)-methylase TrmL, contains SPOUT domain": 'tRNA methylase TrmL',
'Zinc transporter ZupT': 'Zinc transporter ZupT',
'HNH family endonuclease, includes 5-methylcytosine-specific restriction endonuclease McrA': 'HNH endonuclease McrA',
'Cbb3-type cytochrome oxidase, cytochrome c subunit FixO': 'Cytochrome oxidase FixO',
'RecB family exonuclease': 'RecB exonuclease',
'Uncharacterized conserved protein, Alpha-E superfamily': 'Alpha-E conserved protein',
'Glutaredoxin': 'Glutaredoxin',
'Cell division protein CpoB, coordinates peptidoglycan biosynthesis and outer membrane constriction': 'Cell division protein CpoB',
'Predicted lipid-binding transport protein, Tim44 family': 'Tim44 transport protein',
'Uncharacterized conserved protein, contains a C-terminal beta-barrel porin domain': 'Beta-barrel porin protein',
'Sel1-like repeat, TPR-related': 'Sel1 repeat protein',
'Phosphoglycerate dehydrogenase or related dehydrogenase': 'Phosphoglycerate dehydrogenase',
'Zn-dependent oligopeptidase, M3 family': 'Zn oligopeptidase M3',
'ABC-type branched-chain amino acid transport system, permease component': 'ABC BCAA permease',
'Autotransporter adhesin AidA': 'Adhesin AidA',
'Type IV secretion system ATPase VirB4': 'T4SS ATPase VirB4',
'RecA/RadA recombinase': 'RecA recombinase',
'DNA segregation ATPase FtsK/SpoIIIE or related protein': 'DNA segregation ATPase FtsK',
'NADH dehydrogenase, FAD-containing subunit': 'NADH dehydrogenase',
'Outer membrane usher protein FimD/PapC': 'Outer membrane usher FimD/PapC',
'Formylglycine-generating enzyme, required for sulfatase activity, contains SUMF1/FGE domain': 'Formylglycine enzyme',
'Lhr-like helicase': 'Lhr helicase',
'Outer membrane receptor for ferric coprogen and ferric-rhodotorulic acid': 'Ferric coprogen receptor'
}

# Function to apply abbreviation dictionary to row labels
def shorten_label(x):
    for k in abbreviation_dict:
        if k in x:
            return abbreviation_dict[k]
    return x

# Apply function to normalize long labels
rich_factor_df_norm.index = rich_factor_df_norm.index.map(shorten_label)
p_value_df.index = p_value_df.index.map(shorten_label)

# Plot heatmap of normalized rich factors
plt.figure(figsize=(7.6, 11))
ax = sns.heatmap(
    rich_factor_df_norm,
    cmap="Purples",
    linewidths=1,
    linecolor="grey",
    annot=False,
    cbar_kws={"label": "Rich Factor", "fraction": CBAR_WIDTH, "shrink": CBAR_LENGTH}
)

# Overlay significance stars on cells where p-value < 0.05
for i in range(rich_factor_df_norm.shape[0]):
    for j in range(rich_factor_df_norm.shape[1]):
        if p_value_df.iloc[i, j] < 0.05:
            ax.text(j + 0.5, i + 0.95, "*", ha="center", va="center", color="white", fontsize=30, fontweight="bold")

# Set axis labels and tick parameters
ax.set_xlabel("Module", fontsize=AXIS_LABEL_SIZE)
ax.set_ylabel("COG Function", fontsize=AXIS_LABEL_SIZE)
ax.tick_params(axis="x", labelsize=TICK_LABEL_SIZE, rotation=90)
ax.tick_params(axis="y", labelsize=TICK_LABEL_SIZE)

# Customize colorbar appearance
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=CBAR_TICK_SIZE)
cbar.set_label("Scaled Rich Factor", fontsize=CBAR_LABEL_SIZE)

# Save figure as high-resolution PDF
plt.tight_layout()
plt.savefig("Figure5e.pdf", dpi=300, bbox_inches="tight")
plt.show()