## Giotto is a package using python. However, we could install Giotto in R 
# This piple is built up based on https://giottosuite.readthedocs.io/en/latest/subsections/datasets/mouse_visium_kidney.html


# Necessary for installation from R
library(devtools)  # if not installed: install.packages('devtools')
library(remotes)  # if not installed: install.packages('remotes')
remotes::install_github("RubD/Giotto") 

# compilation problems (gfortran)?
# this version does not require C compilation
remotes::install_github("RubD/Giotto@cless") 

if(!"devtools" %in% installed.packages()) {
  install.packages("devtools")
}

devtools::install_github("drieslab/Giotto@suite")

# Ensure Giotto Suite is installed.
if(!"Giotto" %in% installed.packages()) {
  devtools::install_github("drieslab/Giotto@suite")
}

# Ensure GiottoData, a small, helper module for tutorials, is installed.
if(!"GiottoData" %in% installed.packages()) {
  devtools::install_github("drieslab/GiottoData")
}
library(Giotto)
# Ensure the Python environment for Giotto has been installed.
genv_exists = checkGiottoEnvironment()
if(!genv_exists){
  # The following command need only be run once to install the Giotto environment.
  installGiottoEnvironment()
}

#1. Set working directory
results_folder = "C:\\Users\\leh\\OneDrive - TUNI.fi\\Documents\\Data\\capture_area_1_1\\outs\\images_analysis\\Giotto_images"

# Optional: Specify a path to a Python executable within a conda or miniconda
# environment. If set to NULL (default), the Python executable within the previously
# installed Giotto environment will be used.
python_path = "C:\\Users\\leh\\AppData\\Local\\r-miniconda\\envs\\giotto_env\\python.exe" # alternatively, "/local/python/path/python" if desired.

### PART 1: Giotto global instructions and preparations
# creat instruction

## create instructions
install.packages('R.utils')
reticulate::miniconda_uninstall()
installGiottoEnvironment(force_environment = TRUE, force_miniconda = TRUE)
instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = TRUE,
                                  show_plot = FALSE,
                                  python_path = python_path)

## provide path to visium folder
data_path = 'C:\\Users\\leh\\OneDrive - TUNI.fi\\Documents\\Data\\capture_area_1_1\\outs'

### PART 2: Create Giotto object and process data

## directly from visium folder
visium_cap1 = createGiottoVisiumObject(visium_dir = 'C:\\Users\\leh\\OneDrive - TUNI.fi\\Documents\\Data\\capture_area_1_1\\outs',
                                         expr_data = "raw",
                                         png_name = 'tissue_lowres_image.png',
                                         gene_column_index = 2,
                                         instructions = instrs)

## check metadata
pDataDT(visium_cap1)

# check available image names
showGiottoImageNames(visium_cap1) # "image" is the default name

## show aligned image
spatPlot(gobject = visium_cap1, cell_color = 'in_tissue', show_image = T, point_alpha = 0.7)


## subset on spots that were covered by tissue
metadata = pDataDT(visium_cap1)
in_tissue_barcodes = metadata[in_tissue == 1]$cell_ID
visium_cap1 = subsetGiotto(visium_cap1, cell_ids = in_tissue_barcodes)

## filter (this step is necessery if we work directly on raw data file. However, here, filter data is the starting point, therefore, no need for filterGiotto command with threshold 1)
visium_cap1 <- filterGiotto(gobject = visium_cap1,
                              expression_threshold = 1,
                              feat_det_in_min_cells = 30,
                              min_det_feats_per_cell = 600,
                              expression_values = c('raw'),
                              verbose = T)# no work for this step because working in filter

## normalize
visium_cap1 <- normalizeGiotto(gobject = visium_cap1, scalefactor = 6000, verbose = T)

## add gene & cell statistics
visium_cap1 <- addStatistics(gobject = visium_cap1)

## visualize
spatPlot2D(gobject = visium_cap1, show_image = T, point_alpha = 0.7)
spatPlot2D(gobject = visium_cap1, show_image = T, point_alpha = 0.7, cell_color = 'nr_feats', color_as_factor = F)

### PART 3: Dimension reduction
# highly varaiable feature (genes)
visium_cap1 <- calculateHVF(gobject = visium_cap1)

## run PCA on expression values (default)
visium_cap1 <- runPCA(gobject = visium_cap1)
screePlot(visium_cap1, ncp = 30)
plotPCA(gobject = visium_cap1)

## run UMAP and tSNE on PCA space (default)
visium_cap1 <- runUMAP(visium_cap1, dimensions_to_use = 1:30)
plotUMAP(gobject = visium_cap1)
visium_cap1 <- runtSNE(visium_cap1, dimensions_to_use = 1:30)
plotTSNE(gobject = visium_cap1)

### PART4: Cluster
## sNN network (default)
visium_cap1 <- createNearestNetwork(gobject = visium_cap1, dimensions_to_use = 1:30, k = 15)
## Leiden clustering
visium_cap1 <- doLeidenCluster(gobject = visium_cap1, resolution = 0.4, n_iterations = 1000)
plotUMAP(gobject = visium_cap1, cell_color = 'leiden_clus', show_NN_network = T, point_size = 2.5)


### PART5: Co-visualize
# expression and spatial
spatDimPlot(gobject = visium_cap1, cell_color = 'leiden_clus',
            dim_point_size = 2, spat_point_size = 2.5)

spatDimPlot(gobject = visium_cap1, cell_color = 'nr_feats', color_as_factor = F,
            dim_point_size = 2, spat_point_size = 2.5)


### PART 6: cell type marker gene detection 

# Gini markers
gini_markers_subclusters = findMarkers_one_vs_all(gobject = visium_cap1,
                                                  method = 'gini',
                                                  expression_values = 'normalized',
                                                  cluster_column = 'leiden_clus',
                                                  min_featss = 20,
                                                  min_expr_gini_score = 0.5,
                                                  min_det_gini_score = 0.5)
topgenes_gini = gini_markers_subclusters[, head(.SD, 2), by = 'cluster']$feats
topgenes_gini = gini_markers_subclusters[, head(.SD, 2), by = 'cluster']$gene

# violinplot
violinPlot(visium_cap1, feats = unique(topgenes_gini), cluster_column = 'leiden_clus',
           strip_text = 8, strip_position = 'right') # This command is not suitable any more for the new version of Giotto library
violinPlot(visium_cap1, genes = topgenes_gini, cluster_column = 'leiden_clus', strip_text = 8, strip_position = 'right')
violinPlot(visium_cap1, genes = topgenes_gini, cluster_column = 'leiden_clus',
           strip_text = 8, strip_position = 'right',
           save_param = c(save_name = '11-z1-violinplot_gini', base_width = 5, base_height = 10))

# cluster heatmap
# cluster heatmap
plotMetaDataHeatmap(visium_cap1,
                    selected_genes = topgenes_gini,
                    metadata_cols = c('leiden_clus'),
                    x_text_size = 10, y_text_size = 10)

reticulate::miniconda_uninstall()
installGiottoEnvironment(force_environment = TRUE, force_miniconda = TRUE)

# umap plots 
spatFeatPlot2D(visium_cap1,
              expression_values = 'scaled',
              feats = gini_markers_subclusters[, head(.SD, 1), by = 'cluster']$feats,
              cow_n_col = 3, point_size = 1)

# Scran
BiocManager::install("scran")
scran_markers_subclusters = findMarkers_one_vs_all(gobject = visium_cap1,
                                                   method = 'scran',
                                                   expression_values = 'normalized',
                                                   cluster_column = 'leiden_clus')

topgenes_scran = scran_markers_subclusters[, head(.SD, 2), by = 'cluster']$feats

Plot1 = violinPlot(visium_cap1, genes = topgenes_scran,
           cluster_column = 'leiden_clus',
           strip_text = 10, strip_position = 'right')

# cluster heatmap
plotMetaDataHeatmap(visium_cap1, selected_genes = topgenes_scran,
                    metadata_cols = c('leiden_clus'))

# umap plots
dimFeatPlot2D(visium_cap1, expression_values = 'scaled',
              feats = scran_markers_subclusters[, head(.SD, 1), by = 'cluster']$feats,
              cow_n_col = 3, point_size = 1)

### PART 7: cell-type annotation
#Visium spatial transcriptomics does not provide single-cell resolution, making cell type annotation a harder problem. Giotto provides 3 ways to calculate enrichment of specific cell-type signature gene list: PAGE, rank, and hypergeometric test
# cell-type annotation is following Step4_Giotto_analysis which is build based on https://giottosuite.readthedocs.io/en/latest/subsections/datasets/mouse_visium_brain.html

### PART8: spatial grid
visium_cap1 <- createSpatialGrid(gobject = visium_cap1,
                                   sdimx_stepsize = 400,
                                   sdimy_stepsize = 400,
                                   minimum_padding = 0)
spatPlot(visium_cap1, cell_color = 'leiden_clus', show_grid = T,
         grid_color = 'red', spatial_grid_name = 'spatial_grid')

### PART9: Spatial network
## delaunay network: stats + creation
plotStatDelaunayNetwork(gobject = visium_cap1, maximum_distance = 400)

visium_cap1 = createSpatialNetwork(gobject = visium_cap1, minimum_k = 0)
showNetworks(visium_cap1)
spatPlot(gobject = visium_cap1, show_network = T,
         network_color = 'blue', spatial_network_name = 'Delaunay_network')

### PART10: spatial genes
# spatial genes
## kmeans binarization
kmtest = binSpect(visium_cap1)
spatFeatPlot2D(visium_cap1, expression_values = 'scaled',
               feats = kmtest$feats[1:6], cow_n_col = 2, point_size = 1.5)

## rank binarization
ranktest = binSpect(visium_cap1, bin_method = 'rank')
spatFeatPlot2D(visium_cap1, expression_values = 'scaled',
               feats = ranktest$feats[1:6], cow_n_col = 2, point_size = 1.5)

## Spatial co-expression patterns
## spatially correlated genes ##
ext_spatial_genes = kmtest[1:500]$feats

# 1. calculate gene spatial correlation and single-cell correlation
# create spatial correlation object
spat_cor_netw_DT = detectSpatialCorFeats(visium_cap1,
                                         method = 'network',
                                         spatial_network_name = 'Delaunay_network',
                                         subset_feats = ext_spatial_genes) # could not find detectSpatialFeats??? because of version???
spat_cor_netw_DT = detectSpatialPatterns(visium_cap1,
                                         method = 'network',
                                         spatial_network_name = 'Delaunay_network',
                                         subset_feats = ext_spatial_genes) # could not find detectSpatialPatterns??? becaue of version???
# 2. identify most similar spatially correlated genes for one gene
Napsa_top10_genes = showSpatialCorFeats(spat_cor_netw_DT, show_top_feats = 10)
Napsa_top10_genes = showSpatialCorFeats(spat_cor_netw_DT, feats = 'SRSF2', show_top_feats = 10)

spatFeatPlot2D(visium_cap1, expression_values = 'scaled',
               feats = c('ABAT', 'ALDOB', 'TREH', 'FBP1'), point_size = 2)

# 3. cluster correlated genes & visualize
spat_cor_netw_DT = clusterSpatialCorFeats(spat_cor_netw_DT, name = 'spat_netw_clus', k = 8)

heatmSpatialCorFeats(visium_cap1, spatCorObject = spat_cor_netw_DT, use_clus_name = 'spat_netw_clus',
                     save_param = c(save_name = '22-z1-heatmap_correlated_genes', save_format = 'pdf',
                                    base_height = 6, base_width = 8, units = 'cm'),
                     heatmap_legend_param = list(title = NULL))

# 4. rank spatial correlated clusters and show genes for selected clusters
netw_ranks = rankSpatialCorGroups(visium_cap1, spatCorObject = spat_cor_netw_DT, use_clus_name = 'spat_netw_clus',
                                  save_param = c(save_name = '22-z2-rank_correlated_groups',
                                                 base_height = 3, base_width = 5))
# 5. top new spatial cluster
top_netw_spat_cluster = showSpatialCorFeats(spat_cor_netw_DT, use_clus_name = 'spat_netw_clus',
                                            selected_clusters = 6, show_top_feats = 1)

# 5. create metagene enrichment score for clusters
top_netw_spat_cluster = showSpatialCorFeats(spat_cor_netw_DT, use_clus_name = 'spat_netw_clus',
                                            selected_clusters = 6, show_top_feats = 1)
cluster_genes_DT = showSpatialCorFeats(spat_cor_netw_DT, use_clus_name = 'spat_netw_clus', show_top_feats = 1)
cluster_genes = cluster_genes_DT$clus; names(cluster_genes) = cluster_genes_DT$feat_ID

visium_cap1 = createMetafeats(visium_cap1, feat_clusters = cluster_genes, name = 'cluster_metagene')

showGiottoSpatEnrichments(visium_cap1)

spatCellPlot(visium_cap1,
             spat_enr_names = 'cluster_metagene',
             cell_annotation_values = netw_ranks$clusters,
             point_size = 1.5, cow_n_col = 4)

### Part 11: HMRF domains
# HMRF requires a fully connected network!
visium_cap1 = createSpatialNetwork(gobject = visium_cap1, minimum_k = 2, name = 'Delaunay_full')

# spatial genes
my_spatial_genes <- kmtest[1:100]$feats

# do HMRF with different betas
hmrf_folder = paste(results_folder,'/','HMRF_results/')
if(!file.exists(hmrf_folder)) dir.create(hmrf_folder, recursive = T)

# if Rscript is not found, you might have to create a symbolic link, e.g.
# cd /usr/local/bin
# sudo ln -s /Library/Frameworks/R.framework/Resources/Rscript Rscript
HMRF_spatial_genes = doHMRF(gobject = visium_cap1,
                            expression_values = 'scaled',
                            spatial_network_name = 'Delaunay_full',
                            spatial_genes = my_spatial_genes,
                            k = 5,
                            betas = c(0, 1, 6),
                            output_folder = paste0(hmrf_folder, '/', 'Spatial_genes/SG_topgenes_k5_scaled'))

## alternative way to view HMRF results
results = writeHMRFresults(gobject = ST_test,
                           HMRFoutput = HMRF_spatial_genes,
                         k = 5, betas_to_view = seq(0, 25, by = 5))
ST_test = addCellMetadata(ST_test, new_metadata = results, by_column = T, column_cell_ID = 'cell_ID')


## add HMRF of interest to giotto object
visium_cap1 = addHMRF(gobject = visium_cap1,
                        HMRFoutput = HMRF_spatial_genes,
                        k = 5, betas_to_add = c(0, 2),
                        hmrf_name = 'HMRF')

## visualize
spatPlot(gobject = visium_cap1, cell_color = 'HMRF_k5_b.0', point_size = 5)
spatPlot(gobject = visium_cap3, cell_color = 'HMRF_k5_b.2', point_size = 5)

### Export and create Giotto Viewer
# check which annotations are available
combineMetadata(visium_kidney)

# select annotations, reductions and expression values to view in Giotto Viewer
viewer_folder = paste0(results_folder, '/', 'cap3_viewer')

exportGiottoViewer(gobject = visium_cap3,
                   output_directory = viewer_folder,
                   factor_annotations = c('in_tissue',
                                          'leiden_clus'),
                   numeric_annotations = c('nr_feats'),
                   dim_reductions = c('tsne', 'umap'),
                   dim_reduction_names = c('tsne', 'umap'),
                   expression_values = 'scaled',
                   expression_rounding = 2,
                   overwrite_dir = T)
