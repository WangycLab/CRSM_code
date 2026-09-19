# -*- coding: utf-8 -*-
"""
Created on Mon Sep 29 18:53:18 2025
@author: ZHENG XINGHAI
"""

import os, glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "Arial"})

X_LABEL_SIZE, Y_LABEL_SIZE = 30, 30
TICK_LABEL_SIZE, CBAR_WIDTH, CBAR_LENGTH, CBAR_TICK_SIZE, CBAR_LABEL_SIZE = 20, 0.1, 1, 18, 20
crop_order = ["Rice", "Wheat", "Soybean"]

files = glob.glob(os.path.join("COG_enrichment_by_crop", "COG_enrichment_crop_*.csv"))
p_value_dict, rich_factor_dict = {}, {}

for file in files:
    crop = os.path.basename(file).replace("COG_enrichment_crop_", "").replace(".csv", "")
    df = pd.read_csv(file).query("overlap > 0 and rich_factor > 0.01").set_index("COG_ID")
    p_value_dict[crop], rich_factor_dict[crop] = df["p_value"], df["rich_factor"]

p_value_df = pd.DataFrame(p_value_dict).reindex(columns=crop_order).fillna(1)
rich_factor_df = pd.DataFrame(rich_factor_dict).reindex(columns=crop_order).fillna(0)

cog_anno = pd.read_csv("cog-24.def.tab", sep="\t", header=None, usecols=[0, 2], on_bad_lines="skip", engine="python")
cog_anno.columns = ["COG_ID", "Annotation"]
cog_anno_dict = cog_anno.set_index("COG_ID")["Annotation"].to_dict()

rich_factor_df.index = rich_factor_df.index.map(lambda x: cog_anno_dict.get(x, x)).str.strip()
p_value_df.index = p_value_df.index.map(lambda x: cog_anno_dict.get(x, x)).str.strip()

rich_factor_df_norm = rich_factor_df.apply(lambda x: (x - x.min()) / (x.max() - x.min()) if x.max() > x.min() else 0)

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
rich_factor_df_norm, p_value_df = rich_factor_df_norm.iloc[row_order], p_value_df.iloc[row_order]

abbreviation_dict = {
    'Phosphoribosylcarboxyaminoimidazole (NCAIR) mutase': 'NCAIR mutase',
    'Nitric oxide response protein NnrS': 'NO response protein NnrS',
    'Transposase': 'Transposase',
    'Formylglycine-generating enzyme, required for sulfatase activity, contains SUMF1/FGE domain': 'Formylglycine enzyme',
    'Sensor histidine kinase DipB regulating citrate/malate metabolism': 'Histidine kinase DipB',
    'Outer membrane protein TolC': 'Outer membrane TolC',
    'Periplasmic subunit MlaC of the ABC-type intermembrane phospholipid transporter Mla': 'Phospholipid transporter MlaC',
    'Lhr-like helicase': 'Lhr-like helicase',
    'HNH family endonuclease, includes 5-methylcytosine-specific restriction endonuclease McrA': 'HNH endonuclease McrA',
    'Outer membrane receptor for ferric coprogen and ferric-rhodotorulic acid': 'Ferric coprogen receptor',
    'Co-chaperonin GroES (HSP10)': 'Co-chaperonin GroES',
    'Outer membrane usher protein FimD/PapC': 'FimD/PapC usher',
    'ABC-type sugar transport system, ATPase component': 'ABC sugar transporter ATPase',
    'Flagellin and related hook-associated protein FlgL': 'Flagellin FlgL',
    'Adenylosuccinate lyase': 'Adenylosuccinate lyase',
    'Membrane protein insertase Oxa1/YidC/SpoIIIJ': 'Membrane insertase YidC',
    'Periplasmic beta-glucosidase and related glycosidases': 'Periplasmic beta-glucosidase',
    'EntF, seryl-AMP synthase component  of non-ribosomal peptide synthetase': 'NRPS enzyme EntF',
    'Translation elongation factor EF-Tu, a GTPase': 'Translation factor EF-Tu',
    'Multidrug efflux pump subunit AcrA (membrane-fusion protein)': 'Efflux pump AcrA',
    'Acyl-CoA reductase or other NAD-dependent aldehyde dehydrogenase': 'Acyl-CoA reductase',
    'SAM-dependent methyltransferase SmtA (CmoB moved to COG2228)': 'SAM methyltransferase SmtA',
    'Chromosomal replication initiation ATPase DnaA': 'Replication initiator DnaA',
    'Signal transduction histidine kinase': 'Histidine kinase',
    'Outer membrane protein OmpA and related peptidoglycan-associated (lipo)proteins': 'Outer membrane OmpA',
    'Ribosomal protein L34': 'Ribosomal protein L34',
    'ABC-type antimicrobial peptide transport system, permease component': 'ABC antimicrobial permease',
    'Cytoplasmic potassium-binding protein Kbp/XkdP/YgaU, contains LysM domain': 'Kbp potassium-binding protein',
    'Periplasmic deferrochelatase/peroxidase EfeB': 'Periplasmic EfeB'
}

def shorten_label(x):
    return next((v for k, v in abbreviation_dict.items() if k in x), x)

rich_factor_df_norm.index = rich_factor_df_norm.index.map(shorten_label)
p_value_df.index = p_value_df.index.map(shorten_label)

plt.figure(figsize=(1.5, 13))
ax = sns.heatmap(
    rich_factor_df_norm,
    cmap="Purples",
    linewidths=1,
    linecolor="grey",
    annot=False,
    cbar_kws={"label": "Rich Factor", "fraction": CBAR_WIDTH, "shrink": CBAR_LENGTH}
)

for i in range(rich_factor_df_norm.shape[0]):
    for j in range(rich_factor_df_norm.shape[1]):
        if p_value_df.iloc[i, j] < 0.05:
            ax.text(j + 0.5, i + 0.9, "*", ha="center", va="center", color="white", fontsize=40, fontweight="bold")

ax.set_xlabel("Crop", fontsize=X_LABEL_SIZE)
ax.set_ylabel("COG Function", fontsize=Y_LABEL_SIZE)
ax.tick_params(axis="x", labelsize=TICK_LABEL_SIZE, rotation=90)
ax.tick_params(axis="y", labelsize=TICK_LABEL_SIZE)

cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=CBAR_TICK_SIZE)
cbar.set_label("Scaled Rich Factor", fontsize=CBAR_LABEL_SIZE)

plt.tight_layout()
plt.savefig("Figure6d.pdf", dpi=300, bbox_inches="tight")
plt.show()