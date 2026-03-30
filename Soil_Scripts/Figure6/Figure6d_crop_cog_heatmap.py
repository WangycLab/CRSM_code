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

# Ensure fonts remain editable in vector graphics
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'

# Define plot parameters and crop order
AXIS_LABEL_SIZE = 30
TICK_LABEL_SIZE = 20
CBAR_WIDTH = 0.1
CBAR_LENGTH = 1
CBAR_TICK_SIZE = 18
CBAR_LABEL_SIZE = 20
crop_order = ["Rice", "Wheat", "Soybean"]

# Read COG enrichment results for all crops and filter by overlap and rich factor
input_dir = "COG_enrichment_by_crop"
files = glob.glob(os.path.join(input_dir, "COG_enrichment_crop_*.csv"))

p_value_dict = {}
rich_factor_dict = {}
for file in files:
    crop_name = os.path.basename(file).replace("COG_enrichment_crop_", "").replace(".csv", "")
    df = pd.read_csv(file)
    df = df[(df["overlap"] >= 0) & (df["rich_factor"] >= 0.1)]
    p_value_dict[crop_name] = df.set_index("COG_ID")["p_value"]
    rich_factor_dict[crop_name] = df.set_index("COG_ID")["rich_factor"]

# Merge enrichment values into matrices and reorder columns by crop
p_value_df = pd.DataFrame(p_value_dict).fillna(1)
rich_factor_df = pd.DataFrame(rich_factor_dict).fillna(0)
rich_factor_df = rich_factor_df.reindex(columns=crop_order)
p_value_df = p_value_df.reindex(columns=crop_order)

# Read COG annotation file and build dictionary for ID → function mapping
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

# Replace COG IDs in matrices with functional annotations and clean whitespace
rich_factor_df.index = rich_factor_df.index.map(lambda x: cog_anno_dict.get(x, x)).str.strip()
p_value_df.index = p_value_df.index.map(lambda x: cog_anno_dict.get(x, x)).str.strip()

# Normalize rich factor values column-wise to 0–1 for visualization
rich_factor_df_norm = rich_factor_df.copy()
for col in rich_factor_df_norm.columns:
    col_min = rich_factor_df_norm[col].min()
    col_max = rich_factor_df_norm[col].max()
    if col_max > col_min:
        rich_factor_df_norm[col] = (rich_factor_df_norm[col] - col_min) / (col_max - col_min)
    else:
        rich_factor_df_norm[col] = 0

# Cluster rows based on normalized rich factors to determine display order
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

# Define abbreviation dictionary to shorten long COG function names for y-axis
abbreviation_dict = {
    'Phosphoribosylcarboxyaminoimidazole (NCAIR) mutase': 'NCAIR mutase',
    'Cytoplasmic potassium-binding protein Kbp/XkdP/YgaU, contains LysM domain': 'Kbp potassium-binding protein',
    'Periplasmic deferrochelatase/peroxidase EfeB': 'Periplasmic EfeB',
    'Transposase': 'Transposase',
    'Nitric oxide response protein NnrS': 'NO response protein NnrS',
    'Molybdopterin-guanine dinucleotide biosynthesis protein A': 'MGD biosynthesis protein A',
    'NADPH-dependent curcumin reductase CurA': 'CurA reductase',
    'Sensor histidine kinase DipB regulating citrate/malate metabolism': 'Histidine kinase DipB',
    'Ribosomal protein L21': 'Ribosomal protein L21',
    'Formylglycine-generating enzyme, required for sulfatase activity, contains SUMF1/FGE domain': 'Formylglycine enzyme',
    'Outer membrane usher protein FimD/PapC': 'Outer membrane usher FimD/PapC',
    'Lhr-like helicase': 'Lhr-like helicase',
    'HNH family endonuclease, includes 5-methylcytosine-specific restriction endonuclease McrA': 'HNH endonuclease McrA',
    'Adenylosuccinate lyase': 'Adenylosuccinate lyase',
    'Cystathionine beta-lyase/cystathionine gamma-synthase': 'Cystathionine beta/gamma enzyme',
    'Co-chaperonin GroES (HSP10)': 'Co-chaperonin GroES',
    'Outer membrane receptor for ferric coprogen and ferric-rhodotorulic acid': 'Ferric coprogen receptor',
    'ABC-type antimicrobial peptide transport system, permease component': 'ABC antimicrobial permease',
    'Membrane protein insertase Oxa1/YidC/SpoIIIJ': 'Membrane insertase YidC',
    'Periplasmic beta-glucosidase and related glycosidases': 'Periplasmic beta-glucosidase',
    'ABC-type sugar transport system, ATPase component': 'ABC sugar transporter ATPase',
    'Flagellin and related hook-associated protein FlgL': 'Flagellin FlgL',
    'Multidrug efflux pump subunit AcrA (membrane-fusion protein)': 'Efflux pump AcrA',
    'Acyl-CoA reductase or other NAD-dependent aldehyde dehydrogenase': 'Acyl-CoA reductase',
    'Signal transduction histidine kinase': 'Histidine kinase',
    'Outer membrane protein OmpA and related peptidoglycan-associated (lipo)proteins': 'Outer membrane OmpA',
    'SAM-dependent methyltransferase SmtA (CmoB moved to COG2228)': 'SAM methyltransferase SmtA',
    'Chromosomal replication initiation ATPase DnaA': 'Replication initiator DnaA',
    'Periplasmic subunit MlaC of the ABC-type intermembrane phospholipid transporter Mla': 'Phospholipid transporter MlaC',
    'Ribosomal protein L34': 'Ribosomal protein L34',
    'NADP-dependent 3-hydroxy acid dehydrogenase YdfG': '3-hydroxy acid dehydrogenase YdfG',
    'EntF, seryl-AMP synthase component  of non-ribosomal peptide synthetase': 'NRPS enzyme EntF',
    'Outer membrane protein TolC': 'Outer membrane TolC'
}

# Shorten function names using the abbreviation dictionary
def shorten_label(x):
    for k in abbreviation_dict:
        if k in x:
            return abbreviation_dict[k]
    return x

rich_factor_df_norm.index = rich_factor_df_norm.index.map(shorten_label)
p_value_df.index = p_value_df.index.map(shorten_label)

# Plot heatmap of scaled rich factors with significance stars
plt.figure(figsize=(1.8, 13))
ax = sns.heatmap(
    rich_factor_df_norm,
    cmap="Purples",
    linewidths=1,
    linecolor="grey",
    annot=False,
    cbar_kws={"label": "Rich Factor", "fraction": CBAR_WIDTH, "shrink": CBAR_LENGTH}
)

# Add stars for COG functions with p-value < 0.05
for i in range(rich_factor_df_norm.shape[0]):
    for j in range(rich_factor_df_norm.shape[1]):
        if p_value_df.iloc[i, j] < 0.05:
            ax.text(j + 0.5, i + 0.9, "*", ha="center", va="center", color="white", fontsize=40, fontweight="bold")

# Configure axis labels and tick formatting
ax.set_xlabel("Crop", fontsize=AXIS_LABEL_SIZE)
ax.set_ylabel("COG Function", fontsize=AXIS_LABEL_SIZE)
ax.tick_params(axis="x", labelsize=TICK_LABEL_SIZE, rotation=90)
ax.tick_params(axis="y", labelsize=TICK_LABEL_SIZE)

# Customize colorbar appearance
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=CBAR_TICK_SIZE)
cbar.set_label("Scaled Rich Factor", fontsize=CBAR_LABEL_SIZE)

# Save the heatmap figure as a PDF
plt.tight_layout()
plt.savefig("Figure6d.pdf", dpi=300, bbox_inches="tight")
plt.show()