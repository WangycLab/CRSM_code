# CRSM_code

For questions regarding the dataset or analysis, please contact:

- **Yongcheng Wang** - [yongcheng@zju.edu.cn](mailto:yongcheng@zju.edu.cn)
- **Xinghai Zheng** - [xhzheng@zju.edu.cn](mailto:xhzheng@zju.edu.cn)

## Data and Script Preparation

The folder **Soil_Scripts** contains the analysis scripts for **microbial single-cell RNA sequencing (mscRNA-seq) data from crop rhizosphere soils**.

The corresponding **expression matrices** and **species annotation files** for the microbial single-cell RNA sequencing dataset are available at:

https://doi.org/10.6084/m9.figshare.31872934

Please download the dataset and extract it to obtain the folder **Soil_Matrix**.

For subsequent analyses, place the **Soil_Matrix** folder and the **Soil_Scripts** folder in the **same directory path**, as shown below:

```
project_directory/
├── Soil_Scripts
│   ├── Figure3
│   ├── Figure4
│   ├── Figure5
│   └── Figure6
└── Soil_Matrix
```

Once the folders are organized in this structure, the analysis scripts can be executed directly.

---

## Bioinformatics Analysis Workflow

The scripts used to generate each figure are organized into separate folders (`Figure 3`–`Figure 6`) within the **Soil_Scripts** directory.  
The execution order of the scripts for reproducing each figure is described below.

---

## Figure 3

Navigate to the folder:

```
Soil_Scripts/Figure3
```

Run the following script to generate **Figure 3a**:

```
Figure3a_msc_vs_meta.py
```

To generate **Figure 3b–c**, run the scripts in the following order:

```
Figure3bcde_purity_qc.py
Figure3bc_similarity_mds_jsd.py
```

To generate **Figure 3d–e**, run:

```
Figure3de_abundance_profile.py
```

---

## Figure 4

Navigate to the folder:

```
Soil_Scripts/Figure4
```

To generate **Figure 4b**, run:

```
Figure4b_expression_profile.R
Figure4b_expression_similarity.py
```

To generate **Figure 4a, 4c, and 4d**, run the following scripts in sequence:

```
Figure4acd_umap_clustering.R
Figure4a_gene_cell_counts.py
Figure4c_composition_umap.py
Figure4d_sankey_mapping.py
```

These scripts will produce **Figure 4a**, **Figure 4c**, and **Figure 4d**, respectively.

---

## Figure 5

Navigate to the folder:

```
Soil_Scripts/Figure5
```

To generate **Figure 5a**, run the scripts in the following order:

```
Figure5a_cluster_degs.R
Figure5a_cluster_annotation.py
Figure5a_cluster_similarity.py
Figure5a_cluster_feature.py
```

To generate **Figure 5b–c**, run:

```
Figure5bcde_hdWGCNA.R
```

To generate **Figure 5d**, run:

```
Figure5de_module_annotation.py
Figure5d_hubgene_crop_bubble.py
```

To generate **Figure 5e**, run:

```
Figure5e_module_cog_enrichment.py
Figure5e_module_cog_heatmap.py
```

---

## Figure 6

Navigate to the folder:

```
Soil_Scripts/Figure6
```

To generate **Figure 6a**, run:

```
Figure6abcd_single_species_umap.R
```

To generate **Figure 6b**, run:

```
Figure6b_sankey_mapping.py
```

To generate **Figure 6c**, run the following scripts sequentially:

```
Figure6cd_crop_annotation.py
Figure6c_deg_crop_bubble.py
```

To generate **Figure 6d**, run:

```
Figure6d_crop_cog_enrichment.py
Figure6d_crop_cog_heatmap.py
```

The analysis procedures for **Figure 6e–h** are the same as those used for **Figure 6a–d**, except that the analyzed species are replaced. Therefore, the corresponding scripts are not repeated here.
