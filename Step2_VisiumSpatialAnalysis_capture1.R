library(Seurat)
library(SeuratDisk)
### Setup for reference
setwd("C:\\Users\\leh\\OneDrive - TUNI.fi\\Documents\\Data\\cap_3\\outs\\reference_gut")
# Load reference
cortex_sc = readRDS(glue::glue("C:\\Users\\leh\\OneDrive - TUNI.fi\\Documents\\Data\\cap_3\\outs\\reference_gut\\Full_obj_log_counts_soupx_v2_DUOdata.RDS"))

library(dplyr)
library(SummarizedExperiment)
library(SingleCellExperiment)
logcounts(cortex_sc) = assay(cortex_sc, "X")
reference = as.Seurat(cortex_sc, "X")
class(reference)
reference_SCT <- Seurat::SCTransform(reference, assay = NULL) %>%
  Seurat::RunPCA(., verbose = FALSE) %>%
  Seurat::RunUMAP(., dims = 1:30, verbose = FALSE)

### Setup for data
setwd("C:/Users/leh/OneDrive - TUNI.fi/Documents/Data/capture_area_1_1/outs")
data_cap1 = Load10X_Spatial(
  data.dir = "C:/Users/leh/OneDrive - TUNI.fi/Documents/Data/capture_area_1_1/outs",
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
data_cap1_sct = SCTransform(data_cap1, assay = "Spatial", verbose = FALSE)
data_cap1_sct_pca = RunPCA(data_cap1_sct)
data_cap1_sct_pca = RunPCA(data_cap1_sct, assay = "SCT", verbose = FALSE)
data_cap1_sct_pca_cluster = FindNeighbors(data_cap1_sct_pca, reduction = "pca", dims = 1:30)
data_cap1_sct_pca_cluster = FindClusters(data_cap1_sct_pca_cluster, verbose = FALSE)
data_cap1_sct_pca_cluster = RunUMAP(data_cap1_sct_pca_cluster, reduction = "pca", dims = 1:30)
SpatialFeaturePlot(data_cap1_sct_pca_cluster, features = "GBP1", alpha = c(0.1, 1), max.cutoff = "q95")

# Clustering individually
SpatialDimPlot(data_cap1_sct_pca_cluster, cells.highlight = CellsByIdentities(object = data_cap1_sct_pca_cluster, idents = c(0, 1, 2, 3, 4, 5, 6, 7)), facet.highlight = TRUE, ncol = 4)
SpatialDimPlot(data_cap1_sct_pca_cluster, cells.highlight = CellsByIdentities(object = data_cap1_sct_pca_cluster, idents = c(0, 1, 2, 3,4,5,6,7)), facet.highlight = TRUE, ncol = 4)

FeaturePlot( object = data_cap1_sct_pca_cluster, features = c("IFI27", "HLA-C", "APOA1", "APOC3", "DEFA5", "ALPI"), ncol = 3)

### CORTEX 0 and 5
# Subset data from visium spatial data based on cluster 0 and 5 for the right sample in capture_area_1
cortex <- subset(data_cap1_sct_pca_cluster, idents = c(0, 5))
SpatialFeaturePlot(cortex, features = "HMGCS2", alpha = c(0.1, 1), max.cutoff = "q95")
SpatialFeaturePlot(cortex, features = "HMGCS2", alpha = c(0.1, 1), max.cutoff = "q95", pt.size.factor = 1)
# After subsetting, we renormalize cortex
cortex <- SCTransform(cortex, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(verbose = FALSE)
# the annotation is stored in the 'predicted_celltype' column of object metadata
DimPlot(reference_SCT, group.by = "category", label = TRUE)

# match the reference_SCT and cortex for annotating cell types based on clustering
anchors <- FindTransferAnchors(reference = reference_SCT, query = cortex, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = reference_SCT$category, prediction.assay = TRUE,
                                  weight.reduction = cortex[["pca"]], dims = 1:25, k.weight = 60) 

cortex[["predictions"]] <- predictions.assay
DefaultAssay(cortex) <- "predictions"
cortex$predicted.id = GetTransferPredictions(cortex)
Idents(cortex) = "predicted.id"
SpatialFeaturePlot(cortex, features = c("Endothelial", "Epithelial", "Plasma cells", "Mesenchymal", "Neuronal", "Myeloid", "T cells"), pt.size.factor = 1.6, ncol = 3, crop = TRUE, , alpha = c(0.1, 2))
SpatialDimPlot(cortex, cells.highlight = CellsByIdentities(object = cortex, idents = c("Endothelial", "Epithelial", "Mesenchymal", "Myeloid", "Neuronal", "Plasma cells", "T cells", "max")), facet.highlight = TRUE)
SpatialDimPlot(cortex, cells.highlight = CellsByIdentities(object = cortex, idents = c("Endothelial", "Epithelial", "Mesenchymal", "Myeloid", "Neuronal", "Plasma cells", "T cells", "max")), facet.highlight = F, crop = T, image.alpha = 1, label.size = 2, pt.size.factor = 1.6)

### we demonstrate the latter with an implementation of Moran’s I available via FindSpatiallyVariableFeatures() by setting method = 'moransi'
library(ape)
DefaultAssay(cortex) = "SCT"
cortex = FindSpatiallyVariableFeatures(cortex, assay = "SCT", slot = "scale.data", features = VariableFeatures(cortex), selection.method = "moransi", x.cuts = 100, y.cuts = 100)
feature_cortex = HVFInfo(cortex, selection.method = "moransi")

#### sort the data following different columns in feature_cortex for getting the highest features
#1
feature_cortex_subset = feature_cortex[,1]
names(feature_cortex_subset) = rownames(feature_cortex)
feature_cortex_sort = sort(feature_cortex_subset, decreasing = T)
SpatialFeaturePlot(cortex, features = head(names(feature_cortex_sort),
                                              6), ncol = 3, alpha = c(0.1, 1), max.cutoff = "q95")

#2
feature_cortex_subset = feature_cortex[,2]
names(feature_cortex_subset) = rownames(feature_cortex)
feature_cortex_sort = sort(feature_cortex_subset, decreasing = T)
SpatialFeaturePlot(cortex, features = head(names(feature_cortex_sort),
                                           6), ncol = 3, alpha = c(0.1, 1), max.cutoff = "q95")
#3
feature_cortex_subset = feature_cortex[,3]
names(feature_cortex_subset) = rownames(feature_cortex)
feature_cortex_sort = sort(feature_cortex_subset, decreasing = T)
SpatialFeaturePlot(cortex, features = head(names(feature_cortex_sort),
                                           6), ncol = 3, alpha = c(0.1, 1), max.cutoff = "q95")
SpatialFeaturePlot(cortex, features = "HMGCS2", alpha = c(0.1, 1), max.cutoff = "q95")

### CORTEX_of_cluster_0_&_2
cortex <- subset(data_cap1_sct_pca_cluster, idents = c(0,2))
SpatialFeaturePlot(cortex, features = "HMGCS2", alpha = c(0.1, 1), max.cutoff = "q95", pt.size.factor = 1)

### CORTEX_of_cluster_0_&_2_&_5
cortex <- subset(data_cap1_sct_pca_cluster, idents = c(0,2,5))
SpatialFeaturePlot(cortex, features = "HMGCS2", alpha = c(0.1, 1), max.cutoff = "q95", pt.size.factor = 1)

### CORTEX_04_07
# Subset data from visium spatial data based on cluster 4 and 7 for the right sample in capture_area_1
cortex <- subset(data_cap1_sct_pca_cluster, idents = c(4, 7))

# After subsetting, we renormalize cortex
cortex <- SCTransform(cortex, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(verbose = FALSE)
# the annotation is stored in the 'predicted_celltype' column of object metadata
DimPlot(reference_SCT, group.by = "category", label = TRUE)

# match the reference_SCT and cortex for annotating cell types based on clustering
anchors <- FindTransferAnchors(reference = reference_SCT, query = cortex, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = reference_SCT$category, prediction.assay = TRUE,
                                  weight.reduction = cortex[["pca"]], dims = 1:25, k.weight = 60) 

cortex[["predictions"]] <- predictions.assay
DefaultAssay(cortex) <- "predictions"
cortex$predicted.id = GetTransferPredictions(cortex)
Idents(cortex) = "predicted.id"
SpatialFeaturePlot(cortex, features = c("Endothelial", "Epithelial", "Plasma cells", "Mesenchymal", "Neuronal", "Myeloid", "T cells"), pt.size.factor = 1.6, ncol = 3, crop = TRUE, , alpha = c(0.1, 2))
SpatialDimPlot(cortex, cells.highlight = CellsByIdentities(object = cortex, idents = c("Endothelial", "Epithelial", "Mesenchymal", "Myeloid", "Neuronal", "Plasma cells", "T cells", "max")), facet.highlight = TRUE)
SpatialDimPlot(cortex, cells.highlight = CellsByIdentities(object = cortex, idents = c("Endothelial", "Epithelial", "Mesenchymal", "Myeloid", "Neuronal", "Plasma cells", "T cells", "max")), facet.highlight = F, crop = T, image.alpha = 1, label.size = 2, pt.size.factor = 1.6)

### we demonstrate the latter with an implementation of Moran’s I available via FindSpatiallyVariableFeatures() by setting method = 'moransi'
library(ape)
DefaultAssay(cortex) = "SCT"
cortex = FindSpatiallyVariableFeatures(cortex, assay = "SCT", slot = "scale.data", features = VariableFeatures(cortex), selection.method = "moransi", x.cuts = 100, y.cuts = 100)
feature_cortex = HVFInfo(cortex, selection.method = "moransi")

#### sort the data following different columns in feature_cortex for getting the highest features
#1
feature_cortex_subset = feature_cortex[,1]
names(feature_cortex_subset) = rownames(feature_cortex)
feature_cortex_sort = sort(feature_cortex_subset, decreasing = T)
SpatialFeaturePlot(cortex, features = head(names(feature_cortex_sort),
                                           6), ncol = 3, alpha = c(0.1, 1), max.cutoff = "q95")

#2
feature_cortex_subset = feature_cortex[,2]
names(feature_cortex_subset) = rownames(feature_cortex)
feature_cortex_sort = sort(feature_cortex_subset, decreasing = T)
SpatialFeaturePlot(cortex, features = head(names(feature_cortex_sort),
                                           6), ncol = 3, alpha = c(0.1, 1), max.cutoff = "q95")
#3
feature_cortex_subset = feature_cortex[,3]
names(feature_cortex_subset) = rownames(feature_cortex)
feature_cortex_sort = sort(feature_cortex_subset, decreasing = T)
SpatialFeaturePlot(cortex, features = head(names(feature_cortex_sort),
                                           6), ncol = 3, alpha = c(0.1, 1), max.cutoff = "q95")

### Spatial deconvolution using RCTD
# While FindTransferAnchors, we could use spacexr package from GitHub for implementing RCTD
options(timeout = 600000000)
devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)

# Counts, cluster, and spot information is extracted from Seurat query and reference objects to construct Reference and SpatialRNA objects used by RCTD for annotation
library(spacexr)

library(dplyr)
library(SummarizedExperiment)
library(SingleCellExperiment)
cortex_sc = readRDS(glue::glue("C:\\Users\\leh\\OneDrive - TUNI.fi\\Documents\\Data\\cap_3\\outs\\reference_gut\\Full_obj_log_counts_soupx_v2_DUOdata.RDS"))
logcounts(cortex_sc) = assay(cortex_sc, "X")
reference = as.Seurat(cortex_sc, "X")

# set up reference
ref = reference
Idents(ref) <- "celltype"

# extract information to pass to the RCTD Reference function
counts <- ref@meta.data$nCount_originalexp
cluster <- as.factor(ref$category)
names(cluster) <- colnames(ref)
nUMI <- ref$nCount_originalexp
names(nUMI) <- colnames(ref)
reference_new <- Reference(counts, cluster, nUMI)
