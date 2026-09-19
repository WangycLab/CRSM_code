# -*- coding: utf-8 -*-
"""
Created on Fri Jul 25 15:46:16 2025
@author: ZHENG XINGHAI
"""

import os
import glob
import pandas as pd
from scipy.stats import ttest_1samp

# Paths
current_dir = os.getcwd()
parts = current_dir.split(os.sep)
project_root = os.sep.join(parts[:-2]) if parts[-2:] == ["Soil_Scripts", "Figure3"] else current_dir
INPUT_DIR = os.path.join(project_root, "Soil_Matrix")
OUTPUT_DIR = os.path.join(current_dir, "soil_species_filted")
os.makedirs(OUTPUT_DIR, exist_ok=True)

summary_list = []

# Process taxonomy reports
for filepath in glob.glob(os.path.join(INPUT_DIR, "*_sc_taxonomy.report")):
    filename = os.path.basename(filepath)
    base = os.path.splitext(filename)[0]
    parts = base.split("_")
    prefix = "_".join(parts[:2])

    print(f"Processing {prefix}...")

    df = pd.read_csv(filepath, sep="\t")
    n_initial = len(df)

    def test_significant(row):
        _, p = ttest_1samp(
            [row["fraction_total_reads2"], row["fraction_total_reads3"]],
            row["fraction_total_reads"],
            alternative="less"
        )
        return p < 0.05

    df["significant"] = df.apply(test_significant, axis=1)
    df_significant = df.loc[df["significant"]].sort_values("new_est_reads", ascending=False)
    n_selected = len(df_significant)

    output_path = os.path.join(OUTPUT_DIR, f"{prefix}_sc_taxonomy_filted.csv")
    df_significant.to_csv(output_path, index=False)

    print(f"{filename} -> Initial: {n_initial}, Selected: {n_selected}")
    print(f"Saved results to {output_path}")

    summary_list.append({"Sample": prefix, "Initial": n_initial, "Selected": n_selected})

# Summary
summary_df = pd.DataFrame(summary_list)
samples = [f"{t}{batch}" for batch in ["_1", "_2", "_3"] for t in ["Rice", "Soybean", "Wheat"]]
summary_df["Sample"] = pd.Categorical(summary_df["Sample"], categories=samples, ordered=True)
summary_df = summary_df.sort_values("Sample")

summary_csv_path = os.path.join(current_dir, "summary_cell_counts.csv")
summary_df.to_csv(summary_csv_path, index=False)

print(f"Saved summary_cell_counts.csv to {current_dir}")