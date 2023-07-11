library(Seurat)
library(SeuratDisk)
setwd("C:\\Users\\leh\\OneDrive - TUNI.fi\\Documents\\Data\\cap_3\\outs\\reference_gut")
# h5ad format
# step1: convert AnnData object to an h5Seurat file
#Convert("Full_obj_log_counts_soupx_v2.h5ad", dest = "h5seurat", overwrite = TRUE)
Convert("Full_obj_raw_counts_nosoupx_v2.h5ad", dest = "h5seurat", overwrite = TRUE)

# step2: Load h5Seurat file into a Seurat object

seurat_anndata = LoadH5Seurat("Full_obj_raw_counts_nosoupx_v2.h5seurat")

Convert("Full_obj_log_counts_soupx_v2.h5ad", dest = "h5seurat", overwrite = TRUE)
seurat_anndata_log = LoadH5Seurat("Full_obj_log_counts_soupx_v2.h5seurat")
