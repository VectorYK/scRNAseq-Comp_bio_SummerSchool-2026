# scRNAseq-Comp_bio_SummerSchool-2026
Tutorial materials for the Computational Biology Summer School: Single-cell RNA-seq analysis using Seurat.
# Single-cell RNA-seq Tutorial

This repository contains the material for the **Single-cell RNA-seq** practical session conducted as part of the **Computational Biology Summer School**.

## Repository structure

```
.
├── data/
├── scripts/

```

## Requirements

- R (≥ 4.3)
- Seurat

Install the required packages by running

```r
source("scripts/install_packages.R")
```

## Data

Download the tutorial dataset and place the following files inside the `data/` directory.

- `filtered_feature_bc_matrix.h5`
- `metadata.csv`

## Running the tutorial

```r
source("scripts/01_scRNA_pipeline.R")
```

The tutorial covers:

- Data loading
- Quality control
- Normalization
- Feature selection
- PCA
- Clustering
- UMAP visualization
- Marker identification
