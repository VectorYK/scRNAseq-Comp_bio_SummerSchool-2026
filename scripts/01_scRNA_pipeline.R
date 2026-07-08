# ==============================================================================
# Single-Cell RNA-Seq Analysis Tutorial
# Libraries: Seurat
# Input: 'filtered_feature_bc_matrix.h5' and 'metadata.csv'
# ==============================================================================

library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
# ------------------------------------------------------------------------------
# STEP 1: Data Loading and Initial Inspection
# ------------------------------------------------------------------------------

# 1. Load the Count Matrix (.h5)
# Read10X_h5 returns a sparse matrix (genes x cells)
counts_matrix <- Read10X_h5("data/filtered_feature_bc_matrix.h5")

# Check the dimensions (Genes x Cells)
print("Matrix Dimensions:")
print(dim(counts_matrix))

# Peek at the first few barcodes (column names) and genes (row names)
print("First 5 Barcodes:")
print(head(colnames(counts_matrix)))
print("First 5 Genes:")
print(head(rownames(counts_matrix)))

# Optional: Check sparsity (percentage of zeros)
# scRNA-seq data is usually >90% sparse
sparsity <- sum(counts_matrix == 0) / length(counts_matrix)
cat(paste0("Matrix Sparsity: ", round(sparsity * 100, 2), "%\n"))

# 2. Load the Metadata (.csv)
# Set the first column as row names (crucial for matching barcodes)
metadata_df <- read.csv("data/metadata.csv", row.names = 1)

# Check dimensions and preview
print("Metadata Dimensions:")
print(dim(metadata_df))
print(head(metadata_df))

# Check for NA values in the metadata
print(paste("Total NA values in metadata:", sum(is.na(metadata_df))))

# ------------------------------------------------------------------------------
# STEP 2: Create Seurat Object and Quality Control
# ------------------------------------------------------------------------------

# Initialize Seurat object
# This automatically intersects the counts and metadata barcodes
seu_obj <- CreateSeuratObject(
  counts = counts_matrix,
  meta.data = metadata_df,
  project = "MyScRNAAnalysis"
)

# View the object summary
print(seu_obj)

# 1. Calculate Mitochondrial Percentage
# Use pattern = "^MT-" for human or "^mt-" for mouse data
seu_obj[["percent.mt"]] <- PercentageFeatureSet(seu_obj, pattern = "^MT-")

# 2. Visualize QC metrics
# 'percent.mito' assumed to be a column in your custom metadata, 'percent.mt' is calculated above
VlnPlot(seu_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# 3. Apply Filtering (Optional/Commented out for this dataset)
# Adjust these numbers based on the Violin Plots above
# seu_obj <- subset(seu_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
# cat(paste("Cells remaining after QC:", ncol(seu_obj), "\n"))

# ------------------------------------------------------------------------------
# STEP 3: Normalization and Feature Selection
# ------------------------------------------------------------------------------

# 1. Normalize the data (LogNormalize is standard)
seu_obj <- NormalizeData(seu_obj, normalization.method = "LogNormalize", scale.factor = 10000)

# 2. Find Variable Features (Top 2000)
# These genes drive the downstream clustering
seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = 2000)

# Intermediate Output: Top 10 variable genes
top10 <- head(VariableFeatures(seu_obj), 10)
print("Top 10 Variable Features:")
print(top10)

plot1 <- VariableFeaturePlot(seu_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
plot1 + plot2
# ------------------------------------------------------------------------------
# STEP 4: Scaling and PCA
# ------------------------------------------------------------------------------

# 1. Scale data (shifting mean to 0, variance to 1)
# Essential before PCA so high-expression genes don't dominate
all.genes <- rownames(seu_obj)
seu_obj <- ScaleData(seu_obj, features = all.genes)

# 2. Run PCA
seu_obj <- RunPCA(seu_obj, features = VariableFeatures(object = seu_obj))

# Intermediate Output: Inspect PCA loadings
print(seu_obj[["pca"]], dims = 1:5, nfeatures = 5)

# 3. Elbow Plot
# Helps determine how many PCs to use for clustering (usually where the curve flattens)
ElbowPlot(seu_obj)

# 4. Visualize PCA vs Existing Metadata
# Check if cell types are already separating in linear space
DimPlot(seu_obj, reduction = "pca") + NoLegend()
DimPlot(seu_obj, reduction = "pca", group.by = "celltype_major") 

# ------------------------------------------------------------------------------
# STEP 5: Clustering and UMAP
# ------------------------------------------------------------------------------

# 1. Construct the Nearest Neighbor Graph
# We use the first 15 PCs (Adjust based on ElbowPlot)
seu_obj <- FindNeighbors(seu_obj, dims = 1:15)

# 2. Cluster at multiple resolutions
# 0.5 = Coarser clusters, 1.0 = Finer clusters
seu_obj <- FindClusters(seu_obj, resolution = 0.5)
seu_obj <- FindClusters(seu_obj, resolution = 1.0)

# 3. Run UMAP (Non-linear dimensionality reduction)
seu_obj <- RunUMAP(seu_obj, dims = 1:15)

# Intermediate Output: Verify new metadata columns (RNA_snn_res.x)
print(head(seu_obj@meta.data))

# 4. Visualize Clusters vs Metadata

p2 <- DimPlot(seu_obj, reduction = "umap", group.by = "celltype_major", label = TRUE)
p3 <- DimPlot(seu_obj, reduction = "umap", group.by = "celltype_minor")

p2
p3

# ------------------------------------------------------------------------------
# STEP 6: Cell Type Frequency Analysis and Filtering
# ------------------------------------------------------------------------------

# 1. Counts for Cell Type Major
cat("\n--- Counts: Cell Type Major ---\n")
major_counts <- table(seu_obj$celltype_major)
print(major_counts)

# 2. Counts for Cell Type Minor
cat("\n--- Counts: Cell Type Minor ---\n")
minor_counts <- table(seu_obj$celltype_minor)
print(minor_counts)

# 3. Identify Rare Types (threshold < 100 cells)
rare_types <- names(minor_counts[minor_counts < 100])
cat("Types to be relabeled as Undefined:\n")
print(rare_types)

# 4. Create Filtered Metadata Column
# Convert to character to allow string modification
seu_obj$cell_type_minor_filtered <- as.character(seu_obj$celltype_minor)

# Relabel rare types to "Undefined"
seu_obj$cell_type_minor_filtered[seu_obj$cell_type_minor_filtered %in% rare_types] <- "Undefined"

# 5. Verify the changes
cat("\n--- New Counts: cell_type_minor_filtered ---\n")
print(table(seu_obj$cell_type_minor_filtered))

# Visualize with new labels
p1 <- DimPlot(seu_obj, reduction = "umap", group.by = "cell_type_minor_filtered")
p1

# Check a specific gene feature
FeaturePlot(seu_obj, features = "KRT17")

# ------------------------------------------------------------------------------
# STEP 7: Marker Identification
# ------------------------------------------------------------------------------

# Set identity to the new filtered column
Idents(seu_obj) <- "cell_type_minor_filtered"

# 1. Find All Markers
# This compares each cluster against all others
all_markers <- FindAllMarkers(seu_obj, 
                              only.pos = TRUE, 
                              min.pct = 0.25, 
                              logfc.threshold = 0.25)

# Inspect raw results
print(head(all_markers))

# 2. Extract Top 3 Markers per Cluster (Base R)
# Split by cluster
marker_list <- split(all_markers, all_markers$cluster)

# Take top 3 rows from each cluster (assumes sorted by validity)
top3_list <- lapply(marker_list, function(x) head(x, n = 3))

# Recombine into one dataframe
top3_markers <- do.call(rbind, top3_list)

# Print summary
print(top3_markers[, c("cluster", "gene", "avg_log2FC")])

# ------------------------------------------------------------------------------
# STEP 8: Visualization of Markers
# ------------------------------------------------------------------------------

# Get unique genes from top list
genes_to_plot <- unique(top3_markers$gene)

# DotPlot visualization
DotPlot(seu_obj, features = genes_to_plot) + 
  RotatedAxis()
