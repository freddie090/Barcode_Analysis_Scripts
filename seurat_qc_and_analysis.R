################################################################################
# Seurat analysis R script. 
# Freddie Whiting
# Jan 2025
################################################################################

# Description: this script takes command line arguments and subsequently loads
# 10X cell ranger outputs, combines them, filters according to the parameters
# provided. It then normalizes the data, scores for cell-cycling 
# genes, regresses out differences in the difference between cell-cycling 
# scoring genes (s - g2m), and compares a PCA pre- and post- this step.

# Additionally, R images are saved throughout so analysis can be picked up at
# later stages without having to re-run the entire process. All files
# are saved within a directory within the run directory specified with 
# run_regexp.

################################################################################

# Load libraries
library(plyr)
library(dplyr)
library(Seurat)
library(patchwork)
library(purrr)
library(ggplot2)
library(cowplot)
library(limma)
library(MAST)
library(tidyr)
library(msigdbr)
library(fgsea)
library(presto)
library(ComplexHeatmap)
library(RColorBrewer)
library(DoubletFinder)

# # Load in command line arguments: 
args <- commandArgs(trailingOnly = TRUE)
try(if(length(args) != 10) stop("There should be 10 command-line arguments."))

run_regexp <- args[1]
samp_regexp <- args[2]
type_regexp <- args[3]
nfeature_RNA_min <- as.numeric(args[4])
ncount_RNA_min <- as.numeric(args[5])
percent.mt_max <- as.numeric(args[6])
variable.features.n <- as.numeric(args[7])
ncells <- as.numeric(args[8])
n_PCA_dims <- as.numeric(args[9])
cluster_res <- as.numeric(args[10])

################################################################################

# Load in some custom color palettes for plotting. 

# Color palettes for plotting
Hco_pal <- colorRampPalette(brewer.pal("YlGnBu", n = 9)[3:6])(4)
Hdt_pal <- colorRampPalette(brewer.pal("RdPu", n = 9)[6:9])(4)
Sco_pal <- colorRampPalette(brewer.pal("PuBu", n = 9)[3:6])(4)
Sdt_pal <- colorRampPalette(brewer.pal("OrRd", n = 9)[6:9])(4)

Hds_pal <- colorRampPalette(brewer.pal("YlOrBr", n = 9)[3:5])(3)
Sds_pal <- colorRampPalette(brewer.pal("YlOrRd", n = 9)[3:5])(3)

HPOT_col <- colorRampPalette(brewer.pal("Greens", n = 9)[4])(1)
SPOT_col <- colorRampPalette(brewer.pal("Greens", n = 9)[7])(1)

HCT_samp_colvals <- c(
  "POTs" = HPOT_col,
  "CO3_P4_EP1" = Hco_pal[1],
  "CO4_P4_EP1" = Hco_pal[4],
  "DSI_P1_EP1" = Hds_pal[1],
  "DSII_P1_EP1" = Hds_pal[3],
  "DT3_P4_EP1" = Hdt_pal[1],
  "DT4_P4_EP1" = Hdt_pal[4]
)

SW6_samp_colvals <- c(
  "POTs" = SPOT_col,
  "CO1_P4_EP1" = Sco_pal[1],
  "CO3_P4_EP1" = Sco_pal[4],
  "DSI_P1_EP1" = Sds_pal[1],
  "DSII_P1_EP1" = Sds_pal[3],
  "DT1_P4_EP1" = Sdt_pal[1],
  "DT3_P4_EP1" = Sdt_pal[4]
)

# Combined one for combined clustering: 
Comb_samp_colvals <- c(
  "HCTbc_POTs" = HPOT_col,
  "HCTbc_CO3_P4_EP1" = Hco_pal[1],
  "HCTbc_CO4_P4_EP1" = Hco_pal[4],
  "HCTbc_DSI_P1_EP1" = Hds_pal[1],
  "HCTbc_DSII_P1_EP1" = Hds_pal[3],
  "HCTbc_DT3_P4_EP1" = Hdt_pal[1],
  "HCTbc_DT4_P4_EP1" = Hdt_pal[4],
  "SW6bc_POTs" = SPOT_col,
  "SW6bc_CO1_P4_EP1" = Sco_pal[1],
  "SW6bc_CO3_P4_EP1" = Sco_pal[4],
  "SW6bc_DSI_P1_EP1" = Sds_pal[1],
  "SW6bc_DSII_P1_EP1" = Sds_pal[3],
  "SW6bc_DT1_P4_EP1" = Sdt_pal[1],
  "SW6bc_DT3_P4_EP1" = Sdt_pal[4]
)

if(run_regexp == "QR_HCTbc"){
    sample_colvals <- HCT_samp_colvals
} else if (run_regexp == "QR_SW6bc"){
    sample_colvals <- SW6_samp_colvals
} else if (run_regexp == "QR_Comb"){
    sample_colvals <- Comb_samp_colvals
}

################################################################################

# Custom factor orderings to coerce the order plot panels are plotted in. 
# Can use the names from the custom color palettes - already in order. 

# HCTbc sample levels. 
HCT_samp_levels <- names(HCT_samp_colvals)

# HCTbc type levels. 
HCT_type_levels <- unique(sub("(\\D{2}).*", names(HCT_samp_colvals), 
                              replacement = "\\1"))

# SW6bc sample levels. 
SW6_samp_levels <- names(SW6_samp_colvals)

# SW6bc type levels. 
SW6_type_levels <- unique(sub("(\\D{2}).*", names(SW6_samp_colvals), 
                              replacement = "\\1"))

# Comb sample levels. 
Comb_samp_levels <- c(paste0("HCTbc_", HCT_samp_levels),
                      paste0("SW6bc_", SW6_samp_levels))

# Comb type levels. 
Comb_type_levels <- c(paste0("HCTbc_", HCT_type_levels), 
                      paste0("SW6bc_", SW6_type_levels))

# Assign them for this run according to the run_regexp:

if(run_regexp == "QR_HCTbc"){
  sample_levels <- HCT_samp_levels
  type_levels <- HCT_type_levels
} else if(run_regexp == "QR_SW6bc"){
  sample_levels <- SW6_samp_levels
  type_levels <- SW6_type_levels
} else if(run_regexp == "QR_Comb"){
  sample_levels <- Comb_samp_levels
  type_levels <- Comb_type_levels
}

################################################################################

# Go to chosen directory and create Seurat analysis folder
script_dir <- "path/to/script/dir"
cellranger_output_dir <- "path/to/cellranger/output/dir"


# Load the custom functions:
source(paste0(script_dir, "/seurat_qc_and_analysis_functions.R"))

# Go to the cellranger output dir: 
setwd(cellranger_output_dir)

# Use the regexp to go to the correct run directory:
out_dir <- grep(run_regexp, list.files(), value = T)
setwd(out_dir)

if(dir.exists("seurat_analysis") == FALSE) {
  dir.create("seurat_analysis")
  setwd("seurat_analysis")
} else {
  setwd("seurat_analysis")
}

# Merge and filter the Seurat objects: 

seurat_comb_fil <- merge_and_filter_seurat(run_regexp,
                                           samp_regexp,
                                           type_regexp,
                                           ncount_RNA_min,
                                           nfeature_RNA_min,
                                           percent_mt_max) 


# Normalise, scale and perform PCA:

seurat_norm_scal <- normalise_scale_pca(seurat_comb_fil,
                                        run_regexp, 
                                        variable.features.n,
                                        cc.genes)


# Normalization and variance stabilization via SCTransform:

seurat_comb_SCT <- sctransform_workflow(seurat_comb_fil,
                                        run_regexp,
                                        variable.features.n,
                                        ncells,
                                        cc.genes)


# Perform PCA, UMAP and doublet detection and removal: 

seurat_comb_SCT <- pca_umap_doublet_workflow(seurat_comb_SCT,
                                             run_regexp,
                                             n_PCA_dims,
                                             exp_doub=300,
                                             sample_levels=samp_levels,
                                             type_levels=type_levels,
                                             sample_colvals=sample_colvals)


# Run clustering and save QC/summary plots per cluster: 

seurat_comb_SCT <- scRNA_cluster_fun(seurat_comb_SCT,
                                     n_PCA_dims = n_PCA_dims,
                                     cluster_res = cluster_res,
                                     run_regexp,
                                     output_dir = "non_int_analysis_SCT")    


# Calculate proportion of cell per cell-cycle stage for given groups:

scRNA_cell_cycle_prop(seurat_comb_SCT, paste0("SCT_snn_res.", cluster_res))

scRNA_cell_cycle_prop(seurat_comb_SCT, "sample")

scRNA_cell_cycle_prop(seurat_comb_SCT, "type")


# Calculate average expression correlation for given groups: 

scRNA_Avg_Exp_Corr_fun(seurat_comb_SCT, paste0("SCT_snn_res.", cluster_res))

scRNA_Avg_Exp_Corr_fun(seurat_comb_SCT, "sample")

scRNA_Avg_Exp_Corr_fun(seurat_comb_SCT, "type") 



# Calculate differentially expressed genes between all pairwise comparisons or 
# between specific groups: 

scRNA_DE_fun(seurat_comb_SCT, paste0("SCT_snn_res.", cluster_res), 
             all_pairwise = FALSE, output_dir = "non_int_analysis_SCT")

scRNA_DE_fun(seurat_comb_SCT, "sample",
             all_pairwise = FALSE, output_dir = "non_int_analysis_SCT")

scRNA_DE_fun(seurat_comb_SCT, "type",
             all_pairwise = TRUE, output_dir = "non_int_analysis_SCT")


# Perform gene set enrichment analysis using 'fgsea' and 'msigdbr' on either
# all pairwise comparisons or between specific groups:

scRNA_GSE_Analysis_fun(seurat_comb_SCT, paste0("SCT_snn_res.", cluster_res), 
                        all_pairwise = FALSE, output_dir = "non_int_analysis_SCT")  

scRNA_GSE_Analysis_fun(seurat_comb_SCT, "sample",   
                        all_pairwise = FALSE, output_dir = "non_int_analysis_SCT")

scRNA_GSE_Analysis_fun(seurat_comb_SCT, "type", 
                        all_pairwise = FALSE, output_dir = "non_int_analysis_SCT")

scRNA_GSE_Analysis_fun(seurat_comb_SCT, "type", 
                        all_pairwise = TRUE, output_dir = "non_int_analysis_SCT")

################################################################################
