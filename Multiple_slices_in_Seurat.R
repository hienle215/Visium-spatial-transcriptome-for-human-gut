# Working with multiple slices in Seurat
library(Seurat)
data_cap1 = Load10X_Spatial(
  data.dir = "C:/Users/leh/OneDrive - TUNI.fi/Documents/Data/capture_area_1_1/outs",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "slice1",
  filter.matrix = TRUE,
  to.upper = FALSE,
  image = NULL,
)
data_cap2 = Load10X_Spatial(
  data.dir = "C:/Users/leh/OneDrive - TUNI.fi/Documents/Data/capture_area_2_1/outs",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "slice1",
  filter.matrix = TRUE,
  to.upper = FALSE,
  image = NULL,
)
data_cap3 = Load10X_Spatial(
  data.dir = "C:/Users/leh/OneDrive - TUNI.fi/Documents/Data/capture_area_3_1/capture_area_3_1/outs",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "slice1",
  filter.matrix = TRUE,
  to.upper = FALSE,
  image = NULL,
)
data_cap4 = Load10X_Spatial(
  data.dir = "C:/Users/leh/OneDrive - TUNI.fi/Documents/Data/capture_area_4_1/outs",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "slice1",
  filter.matrix = TRUE,
  to.upper = FALSE,
  image = NULL,
)

#In order to work with multiple slices in the same Seurat object, we provide the merge function.
data_cap1 = SCTransform(data_cap1, assay = "Spatial", verbose = FALSE)
data_cap2 = SCTransform(data_cap2, assay = "Spatial", verbose = FALSE)
data_cap3 = SCTransform(data_cap3, assay = "Spatial", verbose = FALSE)
data_cap4 = SCTransform(data_cap4, assay = "Spatial", verbose = FALSE)

data_visium = merge(data_cap1, data_cap2)
data_visium = merge(data_visium, data_cap3)
data_visium = merge(data_visium, data_cap4)

#This then enables joint dimensional reduction and clustering on the underlying RNA expression data.
DefaultAssay(data_visium) <- "SCT"
VariableFeatures(data_visium) <- c(VariableFeatures(data_cap1), VariableFeatures(data_cap2), VariableFeatures(data_cap3), VariableFeatures(data_cap4))
data_visium <- RunPCA(data_visium, verbose = FALSE)
data_visium <- FindNeighbors(data_visium, dims = 1:30)
data_visium <- FindClusters(data_visium, verbose = FALSE)
data_visium <- RunUMAP(data_visium, dims = 1:30)

#Finally, the data can be jointly visualized in a single UMAP plot. SpatialDimPlot() and SpatialFeaturePlot() will by default plot all slices as columns and groupings/features as rows.
DimPlot(data_visium, reduction = "umap", group.by = c("ident", "orig.ident"))
DimPlot(data_visium, reduction = "umap")
SpatialDimPlot(data_visium)
SpatialFeaturePlot(data_visium, features = c("HLA-A", "HMGCS2"))

### Gene expression processing
plot1 <- VlnPlot(data_visium, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot1
plot2 <- SpatialFeaturePlot(data_visium, features = "nCount_Spatial") + theme(legend.position = "right")
plot2
plot1+plot2

#We can then visualize the results of the clustering either in UMAP space (with DimPlot()) or overlaid on the image with SpatialDimPlot().
p1_1 <- DimPlot(data_visium, reduction = "umap", label = TRUE)
p2_1 <- SpatialDimPlot(data_visium, label = TRUE, label.size = 2)
p1_1 + p2_1
SpatialDimPlot(data_visium, cells.highlight = CellsByIdentities(object = data_visium, idents = c(0,1,2,3,4,5)), facet.highlight = T, ncol = 6, label = T, label.size = 2)

# load library
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

cortex <- subset(data_visium, idents = c(4,5,11,8,7))
SpatialFeaturePlot(cortex, features = "HMGCS2", alpha = c(0.1, 1), max.cutoff = "q95", pt.size.factor = 1)

### Annotating cell types

anchors <- FindTransferAnchors(reference = reference_SCT, query = data_visium, normalization.method = "SCT") # 725 anchors
predictions.assay <- TransferData(anchorset = anchors, refdata = reference_SCT$category, prediction.assay = TRUE,
                                  weight.reduction = data_visium[["pca"]], dims = 1:25) 

data_visium[["predictions"]] <- predictions.assay
DefaultAssay(data_visium) <- "predictions"
data_visium$predicted.id = GetTransferPredictions(data_visium)
Idents(data_visium) = "predicted.id"
SpatialFeaturePlot(data_visium, features = c("Endothelial", "Epithelial", "Plasma cells", "Mesenchymal", "Neuronal", "Myeloid", "T cells"), pt.size.factor = 1.6, ncol = 4, crop = TRUE, , alpha = c(0.1, 2))
SpatialFeaturePlot(data_visium, features = c("Endothelial", "Epithelial", "Plasma cells", "Mesenchymal"), pt.size.factor = 1.6, ncol = 4, crop = TRUE, , alpha = c(0.1, 2))
SpatialFeaturePlot(data_visium, features = c("Neuronal", "Myeloid", "T cells"), pt.size.factor = 1.6, ncol = 4, crop = TRUE, , alpha = c(0.1, 2))

SpatialDimPlot(data_visium, cells.highlight = CellsByIdentities(object = data_visium, idents = c("Endothelial", "Epithelial", "Mesenchymal", "Myeloid", "Neuronal", "Plasma cells", "T cells", "max")), facet.highlight = TRUE)
DimPlot(data_visium, reduction = "umap", label = TRUE)

