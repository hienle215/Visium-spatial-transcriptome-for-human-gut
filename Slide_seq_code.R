### SLIDE-SEQ

# Dataset
# Using the SeuratData package for easy data access

setwd("C:/Users/leh/OneDrive - TUNI.fi/Documents/Data/cap_3/outs")
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# Data processing
# The initial preprocessign steps for the bead by gene expression data are similar to other spatial Seurat analyese and to typical scRNA-seq experiment.
# Loading the data
data_cap3 = Load10X_Spatial(
  data.dir = "C:/Users/leh/OneDrive - TUNI.fi/Documents/Data/cap_3/outs",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "slice1",
  filter.matrix = TRUE,
  to.upper = FALSE,
  image = NULL,
)

plot1 <- VlnPlot(data_cap3, features = "nCount_Spatial", pt.size = 0, log = TRUE) + NoLegend()
data_cap3$log_nCount_Spatial <- log(data_cap3$nCount_Spatial)
plot2 <- SpatialFeaturePlot(data_cap3, features = "log_nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

# We then normalize the data using sctransform and perform a standard scRNA-seq dimensionality redcution and clustering workflow

slide.seq = data_cap3
slide.seq <- SCTransform(slide.seq, assay = "Spatial", ncells = 3000, verbose = FALSE)
slide.seq <- RunPCA(slide.seq)
slide.seq <- RunUMAP(slide.seq, dims = 1:30)
slide.seq <- FindNeighbors(slide.seq, dims = 1:30)
slide.seq <- FindClusters(slide.seq, resolution = 0.3, verbose = FALSE)

# We can then visualize the results of the clustering either in UMAP spcace with DimPlot or in the bead coordinate space with SpatialDimPlot function
plot1 <- DimPlot(slide.seq, reduction = "umap", label = TRUE)
plot2 <- SpatialDimPlot(slide.seq, stroke = 0)
plot1 + plot2
SpatialDimPlot(slide.seq, cells.highlight = CellsByIdentities(object = slide.seq, idents = c(0,1,2,3)), facet.highlight = TRUE)

# Visualize QC metrics as a violin plot
VlnPlot(slide.seq, features = c("nFeature_Spatial", "nCount_Spatial"), ncol = 2)

# Featurescater is typically used to visualize feature-feature relationships, but can be used
# For anything calculated by the objet, i.e columns in object metadata, PC scorese etc...
plota= FeatureScatter(slide.seq, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
plota

### Normalized data for annotating gene following scRNA with Seurat object
pbmc = NormalizeData(data_cap3, normalization.method = "LogNormalize")
pbmc = FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
## VST fucntion calculates a variance stabilizing tranformation (VST) from the fitted dispersion-mean relations and then transform the count data (normalized by division by the size factors or normalization factors)
# Identify the 10 most highly variable genes

top10 = head(VariableFeatures(pbmc),10)
top10

# Plot variable features with and without labels
plotb = VariableFeaturePlot(pbmc)
plotb
plotc = LabelPoints(plot = plotb, points = top10, repel = TRUE)
plotc

# all gene
all.genes = rownames(pbmc)
pbmc = ScaleData(pbmc, features = all.genes)
pbmc = RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# Examined and visualize PCA results a feww different ways
print(pbmc[["pca"]], dims = 1:2, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, nfeatures = 15, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(pbmc)

### Find the neighbors and clusters based on identification
pbmc = FindNeighbors(pbmc,dims = 1:30)
pbmc = FindClusters(pbmc, resolution = 0.5)
pbmc = RunUMAP(pbmc, dims = 1:30)
DimPlot(pbmc, reduction = "umap")

# Find all markers for cluster 1
cluster1.markers = FindAllMarkers(pbmc, ident.1 = 2, ident.2 = c(0,1))
VlnPlot(pbmc, features = c)
# Integration with a scRNA-seq reference
# To facilitate cell-type annotation of the slide seq dataset, were are leveraging an existing mouse single cell RNAseq dataset from gut library

library(Seurat)
library(SeuratDisk)
allen_reference <- readRDS("fresh_frozen_RNA_PCA.rds")

# note that setting ncells=3000 normalizes the full dataset but learns noise models on 3k
# cells this speeds up SCTransform dramatically with no loss in performance
library(dplyr)
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)

cortex <- subset(slide.seq, idents = c(1, 2))

# After subsetting, we renormalize cortex
cortex <- SCTransform(cortex, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(verbose = FALSE)
# the annotation is stored in the 'subclass' column of object metadata
DimPlot(allen_reference,  label = TRUE)
# Find anchors
anchors <- FindTransferAnchors(reference = ref, query = slide.seq, normalization.method = "SCT", npcs = 50)
