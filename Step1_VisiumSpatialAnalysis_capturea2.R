
### Overview data
setwd("C:/Users/leh/OneDrive - TUNI.fi/Documents/Data/capture_area_2_1/outs")
library(Seurat)
data_cap2 = Load10X_Spatial(
 data.dir = "C:/Users/leh/OneDrive - TUNI.fi/Documents/Data/capture_area_2_1/outs",
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
plot1 <- VlnPlot(data_cap2, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot1
plot2 <- SpatialFeaturePlot(data_cap2, features = "nCount_Spatial") + theme(legend.position = "right")
plot2
wrap_plots(plot1, plot2)

# SCTransform data
data_cap2_sct = SCTransform(data_cap2, assay = "Spatial", verbose = FALSE)
View(data_cap2_sct)
SpatialFeaturePlot(data_cap2_sct, features = c("IGHG2", "RBP2"))
SpatialFeaturePlot(data_cap2_sct, features = c("HLA-DPA1"))
SpatialFeaturePlot(data_cap2_sct, features = c("GBP1"))
SpatialFeaturePlot(data_cap2_sct, features = "SCG2")

### Gene expression visualization
library(ggplot2)
plot <- SpatialFeaturePlot(data_cap2_sct, features = c("IGHG2", "RBP2")) + theme(legend.text = element_text(size = 0),
                                                               legend.title = element_text(size = 20), legend.key.size = unit(1, "cm"))
plot
jpeg(filename = "C:\\Users\\leh\\OneDrive - TUNI.fi\\Documents\\Data\\capture_area_2_1\\outs\\spatial\\aligned_fiducials.jpg", height = 700, width = 1200, quality = 50)
print(plot)
dev.off()

#The default parameters in Seurat emphasize the visualization of molecular data. However, you can also adjust the size of the spots (and their transparency) to improve the visualization of the histology image, by changing the following parameters:
p1 <- SpatialFeaturePlot(data_cap2_sct, features = c("IGHG2", "RBP2"), pt.size.factor = 1)
p1
p2 <- SpatialFeaturePlot(data_cap2_sct, features = c("IGHG2", "RBP2"), alpha = c(0.1, 1))
p2
p1 + p2

### Dimensional reduction, clustering, and visualization
# We can then proceed to run dimensionality reduction and clustering on the RNA expression data, using the same workflow as we use for scRNA-seq analysis

data_cap2_sct_pca = RunPCA(data_cap2_sct)
data_cap2_sct_pca = RunPCA(data_cap2_sct, assay = "SCT", verbose = FALSE)
data_cap2_sct_pca_cluster = FindNeighbors(data_cap2_sct_pca, reduction = "pca", dims = 1:30)
data_cap2_sct_pca_cluster = FindClusters(data_cap2_sct_pca_cluster, verbose = FALSE)
data_cap2_sct_pca_cluster = RunUMAP(data_cap2_sct_pca_cluster, reduction = "pca", dims = 1:30)

#We can then visualize the results of the clustering either in UMAP space (with DimPlot()) or overlaid on the image with SpatialDimPlot().
p1_1 <- DimPlot(data_cap2_sct_pca_cluster, reduction = "umap", label = TRUE)
p2_1 <- SpatialDimPlot(data_cap2_sct_pca_cluster, label = TRUE, label.size = 3)
p1_1 + p2_1
SpatialDimPlot(data_cap2_sct_pca_cluster, cells.highlight = CellsByIdentities(object = data_cap2_sct_pca_cluster, idents = c(0, 1, 2, 3,4,5,6,7)), facet.highlight = TRUE, ncol = 4, label = TRUE, label.size = 3)

### Interactive plotting: Plot in in the visualization format
SpatialDimPlot(data_cap2_sct_pca_cluster, interactive = TRUE)
SpatialFeaturePlot(data_cap2_sct_pca_cluster, features = "RBP2", interactive = TRUE) # selecting one gene name in features
LinkedDimPlot(data_cap2_sct_pca_cluster)

### Identification of spatially variable features based on cluster 1
de_markers <- FindMarkers(data_cap2_sct_pca_cluster, ident.1 = 1, ident.2 = c(0,2,3,4,5,6,7))
SpatialFeaturePlot(object = data_cap2_sct_pca_cluster, features = rownames(de_markers)[0:6], alpha = c(0.1, 1), ncol = 3)

### Identification of spatially variable features based on cluster 0 and 3
de_markers <- FindMarkers(data_cap2_sct_pca_cluster, ident.1 = c(0,3), ident.2 = c(1,2,4,5,6,7))
SpatialFeaturePlot(object = data_cap2_sct_pca_cluster, features = rownames(de_markers)[0:6], alpha = c(0.1, 1), ncol = 3)
VlnPlot(data_cap2_sct_pca_cluster, features = c("HLA-A", "HLA-DRB1", "TMPRSS15", "HLA-C", "ADA", "APOA1"), ncol = 3)

### Visualization for all clusters
# More details about the identification of variable features
install.packages("Rfast2")
data_cap2_sct_pca_cluster = FindSpatiallyVariableFeatures(data_cap2_sct_pca_cluster, assay = "SCT", features = VariableFeatures(data_cap2_sct_pca_cluster)[1:1000], selection.method = "moransi")

# Now we visualize the expression of the top 6 features identified by this measure (Updated: 15.8.23. the 3 following commands do not work)
top.features = head(SpatiallyVariableFeatures(data_cap2_sct_pca_cluster, selection.method = "moransi"), 3) # error in data frame, no idea where why and how
top.features = head(SpatiallyVariableFeatures(data_cap2_sct_pca_cluster@meta.data, selection.method = "moransi"), 3) # error
top.features = head(SpatiallyVariableFeatures(data_cap2_sct_pca_cluster, selection.method = "markvariogram", decreasing = TRUE), 3)# error

data = data_cap2_sct_pca_cluster

### modifitication solution
top.features = rownames(
  dplyr::slice_min(
    data[["SCT"]]@meta.features,
    moransi.spatially.variable.rank,
    n = 6
  )
) # work well! modified 15.08.2023
SpatialFeaturePlot(data_cap2_sct_pca_cluster, features = top.features, ncol = 3, alpha = c(0.1, 1))
VlnPlot(data_cap2_sct_pca_cluster, features = c("HLA-A", "ADA", "IGHJ6", "HLA-DRB1", "TMPRSS15", "HLA-C"), ncol = 3)


### SUBSET OUT ANATOMICAL REGIONS. Not yet modification 15.08.2023
# As with single-cell objects, you can subset the object to focus on a subset of data. Here, we approximately subset the frontal cortex. This process also facilitates the integration of these data with a cortical scRNA-seq dataset in the next section. First, we take a subset of clusters, and then further segment based on exact positions. After subsetting, we can visualize the cortical cells either on the full image, or a cropped image.
# clusters 1 and 7 will focuse on crypts and villues of left sight samples
cortex <- subset(data_cap2_sct_pca_cluster, idents = c(1, 7))
# now remove additional cells, use SpatialDimPlots to visualize what to remove
# SpatialDimPlot(cortex,cells.highlight = WhichCells(cortex, expression = image_imagerow > 400
# | image_imagecol < 150))
cortex <- subset(cortex, anterior1_imagerow > 300 | anterior1_imagecol < 150, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 275 & anterior1_imagecol > 100, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 250 & anterior1_imagecol > 440, invert = TRUE)

### Slide-seq

plot1 = DimPlot(data_cap2_sct_pca_cluster, reduction = "umap", label = TRUE)
plot2 = SpatialDimPlot(data_cap2_sct_pca_cluster, stroke = 0)

# Clustering individually
SpatialDimPlot(data_cap1_sct_pca_cluster, cells.highlight = CellsByIdentities(object = data_cap3_sct_pca_cluster, idents = c(0, 1, 2, 3)), facet.highlight = TRUE)
FeaturePlot( object = data_cap1_sct_pca_cluster, features = c( "YIF1B", "HLA-DPA1", "MIF-AS1", "SRSF2","MFSD11", "SRSF2", "PCBP1-AS1", "DDTL", "AC009570.2"), ncol = 3)

# Find markers from seurat object
# finding the markers based on clusters
install.packages("BiocManager")
BiocManager::install("multtest")
install.packages("metap")

markers = FindMarkers(object = data_cap3_sct_pca_cluster, ident.1 = 2) # the markers are found here is identitied in cluster 2 and other is compared with remaining cluster. If we change ident.1=1, the gene expression will be found in all cluster 1 and other is compared to
SpatialFeaturePlot(data_cap1_sct_pca_cluster, features = c("DNASE1", "MIF-AS1"))
Idents(FindNeighbors(data_cap1_sct_pca_cluster, dims = 1:10))
Idents(FindClusters(data_cap1_sct_pca_cluster, resolution = 0.5))
str(data_cap1_sct_pca_cluster)
data_cap1$gene
VlnPlot(data_cap1_sct_pca_cluster, features = c("YIF1B", "HLA-DPA1", "MIF-AS1", "SRSF2","MFSD11", "SRSF2", "PCBP1-AS1", "DDTL", "AC009570.2"), ncol = 3)
FeaturePlot(data_cap1_sct_pca_cluster, features = "SRSF2", min.cutoff = "q10")
       
# Extend results to the full datasets

obj <- ProjectData(object = data_cap1_sct_pca_cluster,assay = "RNA",full.reduction = "pca.full",sketched.assay = "sketch",sketched.reduction = "pca",umap.model = "umap",dims = 1:50,refdata = list(cluster_full = "seurat_clusters"))



# Heatmap data of gene expression
DimHeatmap(data_cap2_sct_pca_cluster, dims = 1, cells = 500, balanced = TRUE)

#### NOT YET CHECKED OUT FROM HERE 15.08.2023
# Reman cluster 3 ident
Idents(data_cap1_sct_pca_cluster)
data_cap1_rename = RenameIdents(data_cap1_sct_pca_cluster, "2" = "SRSF2")
DimPlot(data_cap1_rename, reduction = "umap", label = TRUE)

SpatialFeaturePlot(data_cap1_sct, features = "APOB")

### Extend results to the full dataset
obj <- ProjectData(object = data_cap1_sct_pca_cluster,assay = "RNA", full.reduction = "pca.full", sketched.assay = "sketch", sketched.reduction = "pca",umap.model = "umap", dims = 1:50, refdata = list(cluster_full = "seurat_clusters"))
DimPlot(data_cap1_sct_pca_cluster, label = T, label.size = 3) + NoLegend()

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

### Integration with single-cell data

data_ref = Load10X_Spatial(
  data.dir = "C:/Users/leh/OneDrive - TUNI.fi/Documents/Data/cap_3/outs",
  filename = "Full_obj_log_counts_soupx_v2.h5ad",
  assay = "Spatial",
  slice = "slice1",
  filter.matrix = TRUE,
  to.upper = FALSE,
  image = NULL,
)




## Giotto for visium data in R

install.packages('devtools')
library(devtools)
install.packages('remotes')
library(remotes) 
options(buildtools.check = function(action) TRUE ) # to let R environment know about the setup new packages
remotes::install_github("RubD/Giotto") 
# compilation problems (gfortran)?
# this version does not require C compilation
remotes::install_github("RubD/Giotto@cless") 

library(Giotto)

# Set a working directory
#results_folder = '/path/to/directory/'
results_folder = 'C:\\Users\\leh\\OneDrive - TUNI.fi\\Documents\\Data\\cap_3\\outs\\Images_from_R\\Analysis_Giotto'

# set Gioyyo python path
# set python path to your preferred python version path
# set python path to NULL if you want to automatically install (only the 1st time) and use the giotto miniconda environment
python_path = NULL 
if(is.null(python_path)) {
  installGiottoEnvironment()
}

# Part1: Giotto global instruction and preparations
## Create instructions
## create instructions
instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = TRUE,
                                  show_plot = FALSE)

## provide path to visium folder
#data_path = '/path/to/Brain_data/'
data_path = "C:\\Users\\leh\\OneDrive - TUNI.fi\\Documents\\Data\\cap_3\\outs"

# Part2: Create Giotto object and process data
## directly from visium folder
visium_cap3 = createGiottoVisiumObject(visium_dir = data_path, expr_data = "filter", png_name = "tissue_lowres_image.png", gene_column_index = 2, instructions = instrs) # work with filter gene, not with raw file as the tutorial

visium_data = createGiottoVisiumObject(visium_dir = data_path, expr_data = "raw",
                                        png_name = "tissue_lowres_image.png",
                                        gene_column_index = 2, instructions = instrs)

spatPlot(gobject = visium_cap3, cell_color = 'in_tissue', show_image = T, point_alpha = 0.7,
         save_param = list(save_name = '2_a_spatplot_image'))

# Check names
showGiottoImageNames(visium_cap3) # "image" is the default name

# check name
showGiottoImageNames(visium_cap3) # "image" is the default name
# adjust parameters to align image (iterative approach)
visium_cap3 = updateGiottoImage(visium_cap3, image_name = 'image',
                                  xmax_adj = 1300, xmin_adj = 1200,
                                  ymax_adj = 1100, ymin_adj = 1000)

# now it's aligned
spatPlot(gobject = visium_cap3, cell_color = 'in_tissue', show_image = T, point_alpha = 0.7,
         save_param = list(save_name = '2_b_spatplot_image_adjusted'))

## check metadata
pDataDT(visium_cap3)

## compare in tissue with provided jpg
spatPlot(gobject = visium_cap3, cell_color = 'in_tissue', point_size = 2,
         cell_color_code = c('0' = 'lightgrey', '1' = 'blue'),
         save_param = list(save_name = '2_c_in_tissue'))

## subset on spots that were covered by tissue
metadata = pDataDT(visium_cap3)
in_tissue_barcodes = metadata[in_tissue == 1]$cell_ID
visium_cap3 = subsetGiotto(visium_cap3, cell_ids = in_tissue_barcodes)

## filter (Think about it again, do we really need this step if we work in filtered file)
visium_cap3 = filterGiotto(gobject = visium_cap3,
                              expression_threshold = 0.2,
                              gene_det_in_min_cells = 10,
                              min_det_genes_per_cell = 1000,
                              expression_values = c('raw'),
                              verbose = T)

## normalize
visium_cap3 = normalizeGiotto(gobject = visium_cap3, scalefactor = 6000, verbose = T)

## add gene & cell statistics
visium_cap3 = addStatistics(gobject = visium_cap3)

## visualize
spatPlot2D(gobject = visium_cap3, show_image = T, point_alpha = 0.7,
           save_param = list(save_name = '2_d_spatial_locations'))

spatPlot2D(gobject = visium_cap3, show_image = T, point_alpha = 0.7,
           cell_color = 'nr_genes', color_as_factor = F,
           save_param = list(save_name = '2_e_nr_genes'))

### diemnsion reduction
## highly variable genes (HVG)
visium_cap3_dimension <- calculateHVG(gobject = visium_cap3,
                              save_param = list(save_name = '3_a_HVGplot'))

## run PCA on expression values (default)
visium_cap3_dimension <- runPCA(gobject = visium_cap3_dimension, center = TRUE, scale_unit = TRUE)
screePlot(visium_cap3_dimension, ncp = 30, save_param = list(save_name = '3_b_screeplot'))

plotPCA(gobject = visium_cap3_dimension,
        save_param = list(save_name = '3_c_PCA_reduction'))

## run UMAP and tSNE on PCA space (default)
visium_cap3_dimension <- runUMAP(visium_cap3_dimension, dimensions_to_use = 1:10)
plotUMAP(gobject = visium_cap3_dimension,
         save_param = list(save_name = '3_d_UMAP_reduction'))

visium_cap3_dimension <- runtSNE(visium_cap3_dimension, dimensions_to_use = 1:10)
plotTSNE(gobject = visium_cap3_dimension,
         save_param = list(save_name = '3_e_tSNE_reduction'))


## Part4: cluster (can not run it well because the mistakes from tSEN reduction)
## sNN network (default)
visium_cap3_dimension = createNearestNetwork(gobject = visium_cap3_dimension, dimensions_to_use = 1:10, k = 15)
## Leiden clustering
visium_cap3_dimension = doLeidenCluster(gobject = visium_cap3_dimension, resolution = 0.4, n_iterations = 1000)
plotUMAP(gobject = visium_cap3_dimension,
         cell_color = 'leiden_clus', show_NN_network = T, point_size = 2.5,
         save_param = list(save_name = '4_a_UMAP_leiden'))

## Part5: co-visualize

# expression and spatial
spatDimPlot(gobject = visium_cap3_dimension, cell_color = 'leiden_clus',
            dim_point_size = 2, spat_point_size = 2.5,
            save_param = list(save_name = '5_a_covis_leiden'))

spatDimPlot(gobject = visium_cap3_dimension, cell_color = 'nr_genes', color_as_factor = F,
            dim_point_size = 2, spat_point_size = 2.5,
            save_param = list(save_name = '5_b_nr_genes'))

## Part6: cell type marker gene detection
gini_markers_subclusters = findMarkers_one_vs_all(gobject = visium_cap3_dimension,
                                                  method = 'gini',
                                                  expression_values = 'normalized',
                                                  cluster_column = 'leiden_clus',
                                                  min_genes = 20,
                                                  min_expr_gini_score = 0.5,
                                                  min_det_gini_score = 0.5)


### SPOTlight data for deconvolution of spatial transcriptomic data as snRNA-seq cell type profiles to determine spatial interactions

