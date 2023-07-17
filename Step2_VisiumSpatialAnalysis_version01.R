library(Seurat)
library(SeuratDisk)
### Setup for reference
setwd("C:\\Users\\leh\\OneDrive - TUNI.fi\\Documents\\Data\\cap_3\\outs\\reference_gut")
# h5ad format
# step1: convert AnnData object to an h5Seurat file https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185224
#Convert("Full_obj_log_counts_soupx_v2.h5ad", dest = "h5seurat", overwrite = TRUE)
Convert("GSE185224_clustered_annotated_adata_k10_lr0.92_v1.7.h5ad", dest = "h5seurat", overwrite = TRUE)

# step2: Load h5Seurat file into a Seurat object

seurat_anndata = readRDS("all_compartments.RDS")

library(dplyr)
reference_SCT <- reference %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)

### Setup for data
setwd("C:/Users/leh/OneDrive - TUNI.fi/Documents/Data/cap_3/outs")
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
DimPlot(reference_SCT, group.by = "Cluster", label = TRUE)

# match the reference_SCT and cortex for annotating cell types based on clustering
anchors <- FindTransferAnchors(reference = reference_SCT, query = cortex, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = reference_SCT@active.ident, prediction.assay = TRUE,
                                  weight.reduction = cortex[["pca"]], dims = 1:25, k.weight = 18)
cortex[["predictions"]] <- predictions.assay

#Now we get prediction scores for each spot for each class. Of particular interest in the frontal cortex region are the laminar excitatory neurons. Here we can distinguish between distinct sequential layers of these neuronal subtypes, for example:
DefaultAssay(cortex) <- "predictions"
SpatialFeaturePlot(cortex, features = c("Muscularis", "Fibroblasts", "Neural", "Endothelium", "Mesothelium", "Pericytes", "Epithelium", "Immune", "max"), pt.size.factor = 1.6, ncol = 3, crop = TRUE, , alpha = c(0.1, 1))

cortex$predicted.id = GetTransferPredictions(cortex)
Idents(cortex) = "predicted.id"
SpatialDimPlot(cortex, cells.highlight = CellsByIdentities(object = cortex, idents = c("Muscularis", "Fibroblasts", "Neural", "Endothelium", "Mesothelium", "Pericytes", "Epithelium", "Immune", "max")), facet.highlight = TRUE)


# Based on the prediction scores, we can also predict cell types whose location is spatially restricted. 
install.packages("ape", dependencies = T)
library(ape)
DefaultAssay(cortex) = "SCT"
cortex = FindSpatiallyVariableFeatures(cortex, assay = "SCT", slot = "scale.data", features = VariableFeatures(cortex), selection.method = "moransi", x.cuts = 100, y.cuts = 100)

#Now we visualize the expression of the top 6 features identified by Moransi
SpatialFeaturePlot(cortex, features = head(SpatiallyVariableFeatures(cortex, selection.method = "moransi"),6), ncol = 3, alpha = c(0.1,1), max.cutoff = "q95")

### Spatial deconvolution using RCTD
# While FindTransferAnchors, we could use spacexr package from GitHub for implementing RCTD
devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)

# Counts, cluster, and spot information is extracted from Seurat query and reference objects to construct Reference and SpatialRNA objects used by RCTD for annotation
library(spacexr)
