setwd("C:/Users/leh/OneDrive - TUNI.fi/Documents/Data/cap_3/outs")
library(Seurat)
data_cap3 = Load10X_Spatial(
  data.dir = "C:/Users/leh/OneDrive - TUNI.fi/Documents/Data/cap_3/outs",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "slice1",
  filter.matrix = TRUE,
  to.upper = FALSE,
  image = NULL,
)
library(ggplot2)
library(patchwork)
library(dplyr)

reference = readRDS("fresh_frozen_RNA_PCA.rds")

# note that setting ncells=3000 normalizes the full dataset but learns noise models on 3k
# cells this speeds up SCTransform dramatically with no loss in performance
library(dplyr)
reference <- SCTransform(reference, ncells = 3000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)

## Working with visium data
data_cap3_sct = SCTransform(data_cap3, assay = "Spatial", verbose = FALSE)
data_cap3_sct_pca = RunPCA(data_cap3_sct)
data_cap3_sct_pca = RunPCA(data_cap3_sct, assay = "SCT", verbose = FALSE)
data_cap3_sct_pca_cluster = FindNeighbors(data_cap3_sct_pca, reduction = "pca", dims = 1:30)
data_cap3_sct_pca_cluster = FindClusters(data_cap3_sct_pca_cluster, verbose = FALSE)
data_cap3_sct_pca_cluster = RunUMAP(data_cap3_sct_pca_cluster, reduction = "pca", dims = 1:30)

# Clustering individually
SpatialDimPlot(data_cap3_sct_pca_cluster, cells.highlight = CellsByIdentities(object = data_cap3_sct_pca_cluster, idents = c(0, 1, 2, 3)), facet.highlight = TRUE)
FeaturePlot( object = data_cap3_sct_pca_cluster, features = c( "YIF1B", "HLA-DPA1", "MIF-AS1", "SRSF2","MFSD11", "SRSF2", "PCBP1-AS1", "DDTL", "AC009570.2"), ncol = 3)

# Subset data from visium spatial data
cortex <- subset(data_cap3_sct_pca_cluster, idents = c(1, 2))

# After subsetting, we renormalize cortex
cortex <- SCTransform(cortex, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(verbose = FALSE)
# the annotation is stored in the 'predicted_celltype' column of object metadata
DimPlot(reference, group.by = "Predicted_celltype", label = TRUE)

anchors <- FindTransferAnchors(reference = reference, query = cortex, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = reference@meta.data$Predicted_celltype, prediction.assay = TRUE,
                                  weight.reduction = cortex[["pca"]], dims = 1:25, k.weight = 8) # set up k.weight increase, error will be occured

cortex[["predictions"]] <- predictions.assay

#Now we get prediction scores for each spot for each class. Of particular interest in the frontal cortex region are the laminar excitatory neurons. Here we can distinguish between distinct sequential layers of these neuronal subtypes, for example:
DefaultAssay(cortex) <- "predictions"
SpatialFeaturePlot(cortex, features = c("unassigned", "epithelium", "fibroblast", "Bcell", "Tcell", "Glia", "endothelial", "myeloid", "max"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)

# However, cortex is so small subset, maybe it would be better to run directly on data_cap3_pca_cluster

data_cap3_sct_pca_cluster <- SCTransform(data_cap3_sct_pca_cluster, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(verbose = FALSE)
anchors <- FindTransferAnchors(reference = reference, query = data_cap3_sct_pca_cluster, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = reference@meta.data$Predicted_celltype, prediction.assay = TRUE,
                                  weight.reduction = data_cap3_sct_pca_cluster[["pca"]], dims = 1:25, k.weight = 12) # set up k.weight increase, error will be occured
data_cap3_sct_pca_cluster[["predictions"]] <- predictions.assay
DefaultAssay(data_cap3_sct_pca_cluster) <- "predictions"
SpatialFeaturePlot(data_cap3_sct_pca_cluster, features = c("unassigned", "epithelium", "fibroblast", "Bcell", "Tcell", "Glia", "endothelial", "myeloid", "max"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
SpatialFeaturePlot(data_cap3_sct_pca_cluster, features = c("epithelium", "fibroblast", "endothelial", "max"), pt.size.factor = 1.6, ncol = 2, crop = TRUE)

# Based on the prediction scores, we can also predict cell types whose location is spatially restricted. 
install.packages("ape", dependencies = T)
library(ape)
data_cap3_sct_pca_cluster <- FindSpatiallyVariableFeatures(data_cap3_sct_pca_cluster, assay = "predictions", selection.method = "moransi",
                                        features = rownames(data_cap3_sct_pca_cluster), r.metric = 5, slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(data_cap3_sct_pca_cluster, selection.method = "moransi"), 4)
SpatialPlot(object = data_cap3_sct_pca_cluster, features = top.clusters, ncol = 2)


