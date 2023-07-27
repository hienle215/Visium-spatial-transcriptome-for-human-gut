#----------------------------------------------------------------------------Gut atlas survey data------------------------------------------------------------------------------
setwd("C:\\Users\\leh\\OneDrive - TUNI.fi\\Documents\\Data\\cap_3\\outs\\reference_gut")
#sceHI <-  readH5AD("./Full_obj_log_counts_soupx_v2.h5ad", reader = "R") #Memory-heavy, performed once
library(SummarizedExperiment)
BiocManager::install("zellkonverter")
install.packages("C:\\Users\\leh\\OneDrive - TUNI.fi\\Documents\\Data\\cap_3\\outs\\Code_steps\\basilisk_1.11.3.zip", dependencies = TRUE, INSTALL_opts = '--no-lock')
library(basilisk)
BiocManager::install("zellkonverter")
libaray(zellkonverter)
BiocManager::install("HDF5Array")
library(HDF5Array)

sceHI <-  readH5AD("./Full_obj_log_counts_soupx_v2.h5ad", reader = "R") #Memory-heavy, performed once
metadata <- colData(sceHI)
metadataDUO <- metadata[metadata$Region.code == "DUO",]
unique(metadataDUO$category) 
#[1] Endothelial  Epithelial   Plasma cells Mesenchymal  Neuronal     Myeloid      T cells   
unique(metadataDUO$Integrated_05) 
#[1] Mature arterial EC           Paneth                       Microfold cell               TA                           Mature venous EC             Enterocyte                   BEST4+ epithelial           
#[8] Tuft                         IgA plasma cell              LEC1 (ACKR4+)                Contractile pericyte (PLN+)  mLN Stroma (FMO2+)           Stem cells                   LEC6 (ADAMTS4+)             
#[15] LEC5 (CLDN11+)               Adult Glia                   Transitional Stromal 3 (C3+) Macrophages                  gdT                          LYVE1+ Macrophage            Goblet cell                 
#[22] Activated CD4 T              IgG plasma cell              cDC2                         M/X cells (MLN/GHRL+)        TRGV2 gdT           

assays(sceHI)
#List of length 1
#names(1): X

DUO.only <- sceHI[, sceHI$Region.code == "DUO"]
colData(DUO.only)

saveRDS(DUO.only, file = "./Full_obj_log_counts_soupx_v2_DUOdata.rds")
sessionInfo()
