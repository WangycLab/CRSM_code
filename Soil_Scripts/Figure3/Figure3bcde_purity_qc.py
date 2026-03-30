# -*- coding: utf-8 -*-
"""
Created on Fri Jul 25 15:46:16 2025
@author: ZHENG XINGHAI
"""

import os
import glob
import pandas as pd
from scipy.stats import ttest_1samp

# Get current working directory
current_dir = os.getcwd()

# Determine project root
parts = current_dir.split(os.sep)
if parts[-2:] == ["Soil_Scripts", "Figure3"]:
    project_root = os.sep.join(parts[:-2])
else:
    project_root = current_dir

# Define input and output directories relative to project root
INPUT_DIR = os.path.join(project_root, "Soil_Matrix")
OUTPUT_DIR = os.path.join(current_dir, "soil_species_filted")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Initialize list to store summary statistics for all samples
summary_list = []

# Process each taxonomy report file
for filepath in glob.glob(os.path.join(INPUT_DIR, "*_sc_taxonomy.report")):
    filename = os.path.basename(filepath)
    base = os.path.splitext(filename)[0]
    parts = base.split("_")
    prefix = parts[0] + "_" + parts[1]

    print(f"Processing {prefix}...")

    # Read the taxonomy report
    df = pd.read_csv(filepath, sep="\t")
    n_initial = df.shape[0]

    # Perform one-sample t-test to filter significant species
    def test_significant(row):
        other_values = [
            row["fraction_total_reads2"],
            row["fraction_total_reads3"]
        ]
        _, p_value = ttest_1samp(other_values, row["fraction_total_reads"], alternative="less")
        return p_value < 0.05

    df["significant"] = df.apply(test_significant, axis=1)
    df_significant = df[df["significant"]].sort_values(by="new_est_reads", ascending=False)
    n_after_significant = df_significant.shape[0]

    # Save filtered results
    output_path = os.path.join(OUTPUT_DIR, f"{prefix}_sc_taxonomy_filted.csv")
    df_significant.to_csv(output_path, index=False)

    print(f"{filename} -> Initial: {n_initial}, Selected: {n_after_significant}")
    print(f"Saved results to {output_path}")

    # Record summary for this sample
    summary_list.append({
        "Sample": prefix,
        "Initial": n_initial,
        "Selected": n_after_significant
    })

# Combine summaries from all samples into a dataframe
summary_df = pd.DataFrame(summary_list)

# Define sample order for sorting
samples = [f"{t}{batch}" for batch in ["_1", "_2", "_3"] for t in ["Rice", "Soybean", "Wheat"]]
summary_df["Sample"] = pd.Categorical(summary_df["Sample"], categories=samples, ordered=True)
summary_df = summary_df.sort_values("Sample")

# Save the final summary table
summary_csv_path = os.path.join(current_dir, "summary_cell_counts.csv")
summary_df.to_csv(summary_csv_path, index=False)
print(f"Saved summary_cell_counts.csv to {current_dir}")