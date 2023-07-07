### Overview data
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

### Gene expression processing
plot1 <- VlnPlot(data_cap3, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot1
plot2 <- SpatialFeaturePlot(data_cap3, features = "nCount_Spatial") + theme(legend.position = "right")
plot2
wrap_plots(plot1, plot2)

# SCTransform data
data_cap3_sct = SCTransform(data_cap3, assay = "Spatial", verbose = FALSE)
View(data_cap3_sct)
SpatialFeaturePlot(data_cap3_sct, features = c("SRSF2", "MFSD11"))
SpatialFeaturePlot(data_cap3_sct, features = c("HLA-DPA1"))

### Gene expression visualization
library(ggplot2)
plot <- SpatialFeaturePlot(data_cap3_sct, features = c("SRSF2")) + theme(legend.text = element_text(size = 0),
                                                                         legend.title = element_text(size = 20), legend.key.size = unit(1, "cm"))
plot
jpeg(filename = "C:\\Users\\leh\\OneDrive - TUNI.fi\\Documents\\Data\\cap_3\\outs\\spatial\\aligned_fiducials.jpg", height = 700, width = 1200, quality = 50)
print(plot)
dev.off()

#The default parameters in Seurat emphasize the visualization of molecular data. However, you can also adjust the size of the spots (and their transparency) to improve the visualization of the histology image, by changing the following parameters:
p1 <- SpatialFeaturePlot(data_cap3_sct, features = "SRSF2", pt.size.factor = 1)
p1
p2 <- SpatialFeaturePlot(data_cap3_sct, features = "SRSF2", alpha = c(0.1, 1))
p2
p1 + p2

### Dimensional reduction, clustering, and visualization
# We can then proceed to run dimensionality reduction and clustering on the RNA expression data, using the same workflow as we use for scRNA-seq analysis

data_cap3_sct_pca = RunPCA(data_cap3_sct)
data_cap3_sct_pca = RunPCA(data_cap3_sct, assay = "SCT", verbose = FALSE)
data_cap3_sct_pca_cluster = FindNeighbors(data_cap3_sct_pca, reduction = "pca", dims = 1:30)
data_cap3_sct_pca_cluster = FindClusters(data_cap3_sct_pca_cluster, verbose = FALSE)
data_cap3_sct_pca_cluster = RunUMAP(data_cap3_sct_pca_cluster, reduction = "pca", dims = 1:30)

#We can then visualize the results of the clustering either in UMAP space (with DimPlot()) or overlaid on the image with SpatialDimPlot().
p1_1 <- DimPlot(data_cap3_sct_pca_cluster, reduction = "umap", label = TRUE)
p2_1 <- SpatialDimPlot(data_cap3_sct_pca_cluster, label = TRUE, label.size = 3)
p1_1 + p2_1
SpatialDimPlot(data_cap3_sct_pca_cluster, cells.highlight = CellsByIdentities(object = data_cap3_sct_pca_cluster, idents = c(0, 1, 2, 3)), facet.highlight = TRUE, ncol = 3)

### Interactive plotting: Plot in in the visualization format
SpatialDimPlot(data_cap3_sct_pca_cluster, interactive = TRUE)
SpatialFeaturePlot(data_cap3_sct_pca_cluster, features = "SRSF2", interactive = TRUE) # selecting one gene name in features
LinkedDimPlot(data_cap3_sct_pca_cluster)

### Identification of spatially variable features
de_markers <- FindMarkers(data_cap3_sct_pca_cluster, ident.1 = 3, ident.2 = 4)
SpatialFeaturePlot(object = data_cap3_sct_pca_cluster, features = rownames(de_markers)[0:3], alpha = c(0.1, 1), ncol = 3)

# More details about the identification of vavriable features
install.packages("Rfast2")
data_cap3_sct_pca_cluster = FindSpatiallyVariableFeatures(data_cap3_sct_pca_cluster, assay = "SCT", features = VariableFeatures(data_cap3_sct_pca_cluster)[1:1000], selection.method = "moransi")

# Now we visualize the expression of the top 6 features identified by this measure
top.features = head(SpatiallyVariableFeatures(data_cap3_sct_pca_cluster, selection.method = "moransi"), 3) # error in data frame, no idea where why and how
SpatialFeaturePlot(data_cap3_sct_pca_cluster, features = top.features, ncol = 3, alpha = c(0.1, 1))

### Subset out anatomical regions

# As with single-cell objects, you can subset the object to focus on a subset of data. Here, we approximately subset the frontal cortex. This process also facilitates the integration of these data with a cortical scRNA-seq dataset in the next section. First, we take a subset of clusters, and then further segment based on exact positions. After subsetting, we can visualize the cortical cells either on the full image, or a cropped image.

cortex <- subset(data_cap3_sct_pca_cluster, idents = c(1, 2))
# now remove additional cells, use SpatialDimPlots to visualize what to remove
# SpatialDimPlot(cortex,cells.highlight = WhichCells(cortex, expression = image_imagerow > 400
# | image_imagecol < 150))
cortex <- subset(cortex, anterior1_imagerow > 300 | anterior1_imagecol < 150, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 275 & anterior1_imagecol > 100, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 250 & anterior1_imagecol > 440, invert = TRUE)

### Slide-seq

plot1 = DimPlot(data_cap3_sct_pca_cluster, reduction = "umap", label = TRUE)
plot2 = SpatialDimPlot(data_cap3_sct_pca_cluster, stroke = 0)
SpatialDimPlot(data_cap3_sct_pca_cluster, cells.highlight = CellsByIdentities(object = data_cap3_sct_pca_cluster, idents = c(0, 1, 2, 3)), facet.highlight = TRUE)
FeaturePlot( object = data_cap3_sct_pca_cluster, features = c( "YIF1B", "HLA-DPA1", "MIF-AS1", "SRSF2","MFSD11", "SRSF2", "PCBP1-AS1", "DDTL", "AC009570.2"), ncol = 3)

# Find markers from seurat object
# finding the markers based on clusters
install.packages("BiocManager")
BiocManager::install("multtest")
install.packages("metap")

markers = FindMarkers(object = data_cap3_sct_pca_cluster, ident.1 = 2) # the markers are found here is identitied in cluster 2 and other is compared with remaining cluster. If we change ident.1=1, the gene expression will be found in all cluster 1 and other is compared to
SpatialFeaturePlot(data_cap3_sct_pca_cluster, features = c("DNASE1", "MIF-AS1"))
Idents(FindNeighbors(data_cap3_sct_pca_cluster, dims = 1:10))
Idents(FindClusters(data_cap3_sct_pca_cluster, resolution = 0.5))
str(data_cap3_sct_pca_cluster)
data_cap3$gene
VlnPlot(data_cap3_sct_pca_cluster, features = c("YIF1B", "HLA-DPA1", "MIF-AS1", "SRSF2","MFSD11", "SRSF2", "PCBP1-AS1", "DDTL", "AC009570.2"), ncol = 3)
FeaturePlot(data_cap3_sct_pca_cluster, features = "SRSF2", min.cutoff = "q10")

# Extend results to the full datasets

obj <- ProjectData(object = data_cap3_sct_pca_cluster,assay = "RNA",full.reduction = "pca.full",sketched.assay = "sketch",sketched.reduction = "pca",umap.model = "umap",dims = 1:50,refdata = list(cluster_full = "seurat_clusters"))



# Heatmap data of gene expression
DimHeatmap(data_cap3_sct_pca_cluster, dims = 1, cells = 500, balanced = TRUE)

# Reman cluster 3 ident
Idents(data_cap3_sct_pca_cluster)
data_cap3_rename = RenameIdents(data_cap3_sct_pca_cluster, "2" = "SRSF2")
DimPlot(data_cap3_rename, reduction = "umap", label = TRUE)

SpatialFeaturePlot(data_cap3_sct, features = "APOB")

### Extend results to the full dataset
obj <- ProjectData(object = data_cap3_sct_pca_cluster,assay = "RNA", full.reduction = "pca.full", sketched.assay = "sketch", sketched.reduction = "pca",umap.model = "umap", dims = 1:50, refdata = list(cluster_full = "seurat_clusters"))
DimPlot(data_cap3_sct_pca_cluster, label = T, label.size = 3) + NoLegend()

### DE and EnrichR pathway visualization barplot
BiocManager::install("topGO")
levels(data_cap3_sct_pca_cluster)
BiocManager::install("clusterProfiler")
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Hs.eg.db")
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
rownames(data_cap3_sct_pca_cluster)
