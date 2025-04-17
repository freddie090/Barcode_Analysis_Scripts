
#' Merge and Filter Seurat Objects
#'
#' This function reads 10X Genomics data from specified directories, creates 
#' Seurat objects, merges them, and applies filtering based on given thresholds.
#' The resulting filtered Seurat object is saved to an Rdata file, along with 
#' pre- and post-filtering diagnostic plots.
#'
#' @param run_regexp A string used as a prefix for saved files.
#' @param samp_regexp A regular expression to extract sample names from directory names.
#' @param type_regexp A regular expression to extract type names from sample names.
#' @param ncount_RNA_min Minimum threshold for RNA molecule counts per cell.
#' @param nfeature_RNA_min Minimum threshold for the number of genes detected per cell.
#' @param percent_mt_max Maximum threshold for mitochondrial gene expression percentage.
#'
#' @return Saves the filtered Seurat object and diagnostic plots to the working directory.
#' If the Rdata file already exists, it is loaded into the environment.
#' 
#' @examples
#' merge_and_filter_seurat(
#'   run_regexp = "run1",
#'   samp_regexp = "(sample)_.*",
#'   type_regexp = ".*_(type)",
#'   ncount_RNA_min = 500,
#'   nfeature_RNA_min = 200,
#'   percent_mt_max = 5
#' )
merge_and_filter_seurat <- function(
  run_regexp,
  samp_regexp,
  type_regexp,
  ncount_RNA_min,
  nfeature_RNA_min,
  percent_mt_max
) {
  if (!file.exists(paste0(run_regexp, "_seurat_comb_postfil.Rdata"))) {
    setwd("../")
    
    samp_dirs <- grep("seurat_analysis", list.files(), invert = TRUE, value = TRUE)
    samp_names <- sub(samp_regexp, samp_dirs, replacement = "\\1")
    type_names <- sub(type_regexp, samp_names, replacement = "\\1")
    
    list_seurats <- list()
    for (i in seq_along(samp_dirs)) {
      seurat_dat <- Read10X(data.dir = samp_dirs[[i]])
      seurat_obj <- CreateSeuratObject(counts = seurat_dat, min.cells = 5, min.features = 100)
      list_seurats[[i]] <- seurat_obj
      names(list_seurats)[[i]] <- samp_names[[i]]
    }
    
    for (i in seq_along(list_seurats)) {
      list_seurats[[i]] <- RenameCells(list_seurats[[i]], add.cell.id = samp_names[[i]])
      list_seurats[[i]][["percent.mt"]] <- PercentageFeatureSet(list_seurats[[i]], pattern = "^MT-")
      list_seurats[[i]][["sample"]] <- samp_names[[i]]
      list_seurats[[i]][["type"]] <- type_names[[i]]
    }
    
    seurat_comb <- purrr::reduce(list_seurats, merge)
    max_ncount <- max(seurat_comb@meta.data$nCount_RNA)
    max_nfeat <- max(seurat_comb@meta.data$nFeature_RNA)
    rm(list_seurats)
    
    prefil_plot <- ggplot(data = seurat_comb@meta.data, 
                          aes(x = nCount_RNA, y = nFeature_RNA, colour = percent.mt)) +
      geom_point() + stat_smooth(method = lm) +
      geom_vline(xintercept = ncount_RNA_min, linetype = "dashed") +
      geom_hline(yintercept = nfeature_RNA_min, linetype = "dashed") +
      scale_x_log10(name = "Molecules/Cell", limits = c(100, max_ncount + 1000)) +
      scale_y_log10(name = "Genes/Cell", limits = c(10, max_nfeat + 1000)) +
      scale_colour_viridis_c(name = "% Mito.") + theme_minimal() +
      facet_wrap(~sample, ncol = 4) + ggtitle("Pre-Filtering")
    
    seurat_comb_fil <- subset(seurat_comb, 
                              subset = (nFeature_RNA > nfeature_RNA_min) & 
                                (nCount_RNA > ncount_RNA_min) & 
                                (percent.mt < percent_mt_max))
    
    postfil_plot <- ggplot(data = seurat_comb_fil@meta.data, 
                           aes(x = nCount_RNA, y = nFeature_RNA, colour = percent.mt)) +
      geom_point() + stat_smooth(method = lm) +
      geom_vline(xintercept = ncount_RNA_min, linetype = "dashed") +
      geom_hline(yintercept = nfeature_RNA_min, linetype = "dashed") +
      scale_x_log10(name = "Molecules/Cell", limits = c(100, max_ncount + 1000)) +
      scale_y_log10(name = "Genes/Cell", limits = c(10, max_nfeat + 1000)) +
      scale_colour_viridis_c(name = "% Mito.") + theme_minimal() +
      facet_wrap(~sample, ncol = 4) + ggtitle("Post-Filtering")
    
    setwd("seurat_analysis")
    n_samps <- length(unique(seurat_comb_fil@meta.data$sample))
    ggsave2(paste0(run_regexp, "_filter_prefil_plot.pdf"), plot = prefil_plot, 
            height = 5 * ceiling(n_samps / 4), width = 10, bg = "white")
    ggsave2(paste0(run_regexp, "_filter_postfil_plot.pdf"), plot = postfil_plot, 
            height = 5 * ceiling(n_samps / 4), width = 10, bg = "white")
    
    save(seurat_comb_fil, file = paste0(run_regexp, "_seurat_comb_postfil.Rdata"))
  } else {
    load(file = paste0(run_regexp, "_seurat_comb_postfil.Rdata"))
  }
}


#' Normalise, Scale, and Perform PCA on Seurat Object
#'
#' This function normalises data, scores for cell cycling, identifies highly variable genes,
#' scales the data, and performs PCA. It saves the processed Seurat object and PCA plots to files.
#'
#' @param seurat_comb_fil A Seurat object to be normalised and processed.
#' @param run_regexp A string used as a prefix for saved files.
#' @param variable.features.n Number of variable features to identify.
#' @param cc.genes A list containing cell cycle gene sets.
#'
#' @return Saves the processed Seurat object and PCA plot. If the Rdata file already exists, it is loaded into the environment.
#'
#' @examples
#' normalise_scale_pca(
#'   seurat_comb_fil = seurat_object,
#'   run_regexp = "run1",
#'   variable.features.n = 2000,
#'   cc.genes = cc.genes
#' )
normalise_scale_pca <- function(seurat_comb_fil, run_regexp, variable.features.n, cc.genes) {
  if (!file.exists(paste0(run_regexp, "_seurat_comb_norm_scal.Rdata"))) {
    
    # Normalise the data
    seurat_norm <- NormalizeData(seurat_comb_fil, normalization.method = "LogNormalize",
                                 scale.factor = 10000)
    
    # Score for Cell Cycling
    seurat_norm <- CellCycleScoring(seurat_norm, 
                                    g2m.features = cc.genes$g2m.genes,
                                    s.features = cc.genes$s.genes,
                                    set.ident = TRUE)
    
    # Identify most highly variable genes
    seurat_norm <- FindVariableFeatures(seurat_norm,
                                        selection.method = "vst",
                                        nfeatures = variable.features.n)
    
    # Scale the counts
    seurat_norm_scal <- ScaleData(seurat_norm)
    
    # Perform the initial PCA
    seurat_norm_scal <- RunPCA(seurat_norm_scal)
    
    # Plot PCA, split by cell cycle phase
    PCA_1 <- DimPlot(seurat_norm_scal,
                     reduction = "pca",
                     group.by = "sample",
                     split.by = "Phase")
    
    # Save the PCA plot
    ggsave2("PCA_pre_cc_reg.pdf", plot = PCA_1, height = 5, width = 10)
    
    # Save the processed Seurat object
    save(seurat_norm_scal, file = paste0(run_regexp, "_seurat_comb_norm_scal.Rdata"))
  } else {
    # Load existing file
    load(file = paste0(run_regexp, "_seurat_comb_norm_scal.Rdata"))
  }
}


#' Perform SCTransform Workflow on a Seurat Object
#'
#' This function applies the SCTransform workflow to a Seurat object. It scores 
#' for cell cycling, calculates cell cycle differences, and applies SCTransform 
#' while regressing out the cell cycle difference. The processed Seurat object 
#' is saved to an Rdata file.
#'
#' @param seurat_comb_fil A Seurat object to be processed.
#' @param run_regexp A string used as a prefix for saved files.
#' @param variable.features.n Number of variable features to identify.
#' @param ncells Number of cells to use for the SCTransform step.
#' @param cc.genes A list containing cell cycle gene sets.
#'
#' @return Saves the processed Seurat object to an Rdata file. If the file already exists, it is loaded into the environment.
#'
#' @examples
#' sctransform_workflow(
#'   seurat_comb_fil = seurat_object,
#'   run_regexp = "run1",
#'   variable.features.n = 3000,
#'   ncells = 5000,
#'   cc.genes = cc.genes
#' )
sctransform_workflow <- function(seurat_comb_fil, run_regexp, variable.features.n, ncells, cc.genes) {
  if (!file.exists(paste0(run_regexp, "_seurat_SCT_comb.Rdata"))) {
    
    # Create a copy of the input Seurat object
    seurat_comb_SCT <- seurat_comb_fil
    
    # Score for cell cycling
    seurat_comb_SCT <- CellCycleScoring(seurat_comb_SCT, 
                                        g2m.features = cc.genes$g2m.genes,
                                        s.features = cc.genes$s.genes)
    
    # Calculate cell cycle difference
    seurat_comb_SCT$CC.Difference <- seurat_comb_SCT$S.Score - seurat_comb_SCT$G2M.Score
    
    # Run SCTransform
    seurat_comb_SCT <- SCTransform(seurat_comb_SCT,
                                   vars.to.regress = "CC.Difference",
                                   variable.features.n = variable.features.n,
                                   ncells = ncells)
    
    # Save the processed Seurat object
    save(seurat_comb_SCT, file = paste0(run_regexp, "_seurat_SCT_comb.Rdata"))
  } else {
    # Load existing file
    load(file = paste0(run_regexp, "_seurat_SCT_comb.Rdata"))
  }
}


#' PCA, UMAP, and Doublet Detection Workflow for Seurat Objects
#'
#' This function performs PCA, UMAP, doublet detection and removal on a Seurat object. 
#' It generates relevant plots and saves the updated Seurat object and plots to files.
#'
#' @param seurat_comb_SCT A Seurat object after SCTransform.
#' @param run_regexp A string used as a prefix for saved files.
#' @param n_PCA_dims Number of PCA dimensions to use for UMAP and doublet detection (defaults to 40).
#' @param exp_doub Expected number of doublets for DoubletFinder (defaults to 300).
#' @param sample_levels A vector specifying custom factor orderings for samples (optional).
#' @param type_levels A vector specifying custom factor orderings for types (optional).
#' @param sample_colvals A named list specifying custom colour palettes for samples (optional).
#'
#' @return Saves the updated Seurat object and various plots. Returns the processed Seurat object.
#'
#' @examples
#' pca_umap_doublet_workflow(
#'   seurat_comb_SCT = seurat_object,
#'   run_regexp = "QR_Comb",
#'   n_PCA_dims = 40,
#'   exp_doub = 300,
#'   sample_levels = sample_levels,
#'   type_levels = type_levels,
#'   sample_colvals = sample_colvals
#' )
pca_umap_doublet_workflow <- function(seurat_comb_SCT, run_regexp, n_PCA_dims = 40, exp_doub = 300,
                                      sample_levels = NULL, type_levels = NULL,
                                      sample_colvals = NULL) {
  # Run PCA
  seurat_comb_SCT <- RunPCA(object = seurat_comb_SCT)
  
  # Run UMAP
  seurat_comb_SCT <- RunUMAP(seurat_comb_SCT, dims = 1:n_PCA_dims, reduction = "pca")
  
  # PCA plot
  PCA_SCT <- DimPlot(seurat_comb_SCT, reduction = "pca", group.by = "sample", split.by = "Phase")
  
  # Doublet detection
  if (!"doublet_finder" %in% colnames(seurat_comb_SCT@meta.data)) {
    seurat_comb_SCT_split <- SplitObject(seurat_comb_SCT, split.by = "sample")
    
    # Doublet detection function
    run_doublet_detec <- function(seurat_split) {
      for (i in seq_along(seurat_split)) {
        seurat_samp <- seurat_split[[i]]
        sweep.list <- paramSweep_v3(seurat_samp, PCs = 1:n_PCA_dims, sct = TRUE)
        sweep.stats <- summarizeSweep(sweep.list)
        bcmvn <- find.pK(sweep.stats)
        optimal.pk <- as.numeric(levels(bcmvn$pK))[which.max(bcmvn$BCmetric)]
        
        seurat_samp <- doubletFinder_v3(seu = seurat_samp, PCs = 1:n_PCA_dims, 
                                        pK = optimal.pk, nExp = exp_doub, sct = TRUE)
        colnames(seurat_samp@meta.data)[grep("DF.classifications", colnames(seurat_samp@meta.data))] <- "doublet_finder"
        seurat_split[[i]] <- seurat_samp
      }
      return(seurat_split)
    }
    
    seurat_comb_SCT_split <- run_doublet_detec(seurat_comb_SCT_split)
    seurat_comb_SCT <- purrr::reduce(seurat_comb_SCT_split, merge)
  }
  
  # Remove doublets
  seurat_comb_SCT <- subset(seurat_comb_SCT, subset = (doublet_finder == "Singlet"))
  
  # Apply custom sample/type ordering if provided
  if (!is.null(sample_levels)) {
    seurat_comb_SCT$sample <- factor(seurat_comb_SCT$sample, levels = samp_levels)
    seurat_comb_SCT$type <- factor(seurat_comb_SCT$type, levels = type_levels)
  }
  
  # Generate UMAP plots
  UMAP_plot1 <- DimPlot(seurat_comb_SCT, group.by = "sample", pt.size = 1.0)
  UMAP_plot2 <- DimPlot(seurat_comb_SCT, group.by = "sample", split.by = "type", ncol = 4, pt.size = 1.0)
  
  # Apply custom colour palettes if provided
  if (!is.null(sample_colvals)) {
    UMAP_plot1 <- UMAP_plot1 + scale_colour_manual(values = alpha(sample_colvals, 0.6))
    UMAP_plot2 <- UMAP_plot2 + scale_colour_manual(values = alpha(sample_colvals, 0.6))
  }
  
  # Save plots
  n_type <- length(unique(seurat_comb_SCT$type))
  ggsave2("PCA_SCT_Assay.pdf", plot = PCA_SCT, height = 5, width = 7)
  ggsave2("UMAP1_SCT_Assay.pdf", plot = UMAP_plot1, height = 5, width = 7)
  ggsave2("UMAP2_SCT_Assay.pdf", plot = UMAP_plot2, 
          height = 5 * ceiling(n_type / 4), 
          width = 22 - (22 * (n_type < 4)) + ((n_type / 4) * 22) * (n_type < 4))
  
  # Save the updated Seurat object
  save(seurat_comb_SCT, file = paste0(run_regexp, "_seurat_SCT_comb.Rdata"))
  
  return(seurat_comb_SCT)
}


#' Clustering Workflow for Seurat Objects
#'
#' This function performs clustering, generates various QC plots, and evaluates potential
#' drivers of clustering on a Seurat object. It supports SCTransformed, non-integrated data.
#'
#' @param seurat_object A Seurat object to be clustered.
#' @param n_PCA_dims Number of PCA dimensions to use for clustering.
#' @param cluster_res Base resolution for clustering.
#' @param run_regexp A string used as a prefix for saved files and to identify datasets.
#' @param output_dir Directory to save the clustering outputs (default: "non_int_analysis_SCT").
#'
#' @return Saves clustering plots and data, and returns the clustered Seurat object.
#'
#' @examples
#' clustered_seurat <- scRNA_cluster_fun(
#'   seurat_object = seurat_comb_SCT,
#'   n_PCA_dims = 20,
#'   cluster_res = 0.5,
#'   run_regexp = "run1"
#' )
scRNA_cluster_fun <- function(seurat_object, n_PCA_dims = 20, cluster_res = 0.5,
                              run_regexp, output_dir = "non_int_analysis_SCT") {
  
  # Create or move to the output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  setwd(output_dir)
  
  # Identify number of samples and types
  n_samps <- length(unique(seurat_object@meta.data$sample))
  n_type <- length(unique(seurat_object@meta.data$type))
  
  # Build the K-nearest neighbour graph
  seurat_object <- FindNeighbors(object = seurat_object, dims = 1:n_PCA_dims)
  
  # Clustering resolutions
  cluster_res_range <- seq(0.2, cluster_res + 0.4, 0.2)
  
  # Perform clustering at various resolutions
  seurat_object <- FindClusters(object = seurat_object, resolution = cluster_res_range)
  
  # Generate clustering comparison plots
  clust_range_names <- paste0("SCT_snn_res.", cluster_res_range)
  clustering_comp_plots <- lapply(clust_range_names, function(res) {
    Idents(seurat_object) <- res
    DimPlot(seurat_object, reduction = "umap", label = TRUE, label.size = 6) + ggtitle(res)
  })
  
  clustering_comp_plot <- cowplot::plot_grid(plotlist = clustering_comp_plots)
  ggsave2("Clustering_Res_Comp_UMAPS.pdf", plot = clustering_comp_plot, height = 12, width = 20)
  
  # Set chosen resolution for downstream analysis
  Idents(seurat_object) <- paste0("SCT_snn_res.", cluster_res)
  
  # QC: Number of cells per cluster
  n_cells <- FetchData(seurat_object, vars = c("ident", "sample", "type")) %>%
    dplyr::count(ident, sample, type) %>%
    tidyr::spread(ident, n)
  
  n_cells_df <- n_cells %>%
    tidyr::gather(key = "cluster", value = "count", -sample, -type) %>%
    replace(is.na(.), 0) %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(prop = count / sum(count)) %>%
    data.frame()
  
  n_cells_plot <- ggplot(n_cells_df, aes(x = sample, y = prop, fill = as.factor(cluster))) +
    geom_col(colour = "black", alpha = 0.6, width = 0.95) +
    theme_minimal() +
    xlab("\nSample") +
    ylab("Proportion\n") +
    theme(text = element_text(size = 22),
          panel.grid = element_blank(),
          axis.ticks.x = element_line(),
          axis.text.x = element_text(angle = 90, hjust = 1))
  
  # Apply custom colour palette if applicable
  n_groups <- length(unique(n_cells_df$cluster))
  if (n_groups <= 12) {
    col_vec <- sample(colorRampPalette(brewer.pal("Paired", n = n_groups))(n_groups), n_groups)
    names(col_vec) <- 0:(n_groups - 1)
    n_cells_plot <- n_cells_plot + scale_fill_manual(values = col_vec)
  }
  
  ggsave2("Cells_per_cluster.pdf", plot = n_cells_plot, height = 8, width = 9, bg = "white")
  
  # Generate UMAP plots split by sample and type
  UMAP_by_samp <- DimPlot(seurat_object, reduction = "umap", label = TRUE, label.size = 6,
                          split.by = "sample", ncol = n_samps, pt.size = 2)
  UMAP_by_type <- DimPlot(seurat_object, reduction = "umap", label = TRUE, label.size = 6,
                          split.by = "type", ncol = 4, pt.size = 2)
  
  if (n_groups <= 12) {
    UMAP_by_samp <- UMAP_by_samp + scale_colour_manual(values = alpha(col_vec, 0.6)) + theme_void()
    UMAP_by_type <- UMAP_by_type + scale_colour_manual(values = alpha(col_vec, 0.6))
  }
  
  ggsave2("Clustering_UMAPS_by_sample.pdf", plot = UMAP_by_samp, height = 5, width = 5.0 * n_samps)
  ggsave2("Clustering_UMAPS_by_type.pdf", plot = UMAP_by_type, 
          height = 5 * (ceiling(n_type / 4)), 
          width = 21 - (21 * (n_type < 4)) + ((n_type / 4) * 21) * (n_type < 4))
  
  # Save the clustered data
  save(seurat_object, file = paste0(run_regexp, "_seurat_SCT_clust_comb.Rdata"))
  
  setwd("../")
  return(seurat_object)
}


#' Proportion of Cells per Cell-Cycle Stage
#'
#' This function calculates and visualises the proportion of cells in each cell-cycle stage 
#' (G1, S, and G2/M) for a specified grouping in a Seurat object. The results are saved as a plot.
#'
#' @param seurat_object A Seurat object containing cell-cycle stage information in the "Phase" column.
#' @param group_diff The grouping variable (e.g., "sample", "type") to calculate proportions by.
#' @param output_dir Directory to save the output plot (default: "non_int_analysis_SCT").
#'
#' @return Saves a bar plot showing the proportion of cells in each cell-cycle stage grouped by the specified variable.
#'
#' @examples
#' scRNA_cell_cycle_prop(
#'   seurat_object = seurat_comb_SCT,
#'   group_diff = "sample"
#' )
scRNA_cell_cycle_prop <- function(seurat_object, group_diff, output_dir = "non_int_analysis_SCT") {
  
  # Create or move to the output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  setwd(output_dir)
  
  # Set group identity
  Idents(seurat_object) <- group_diff
  seurat_object@meta.data$group <- Idents(seurat_object)
  
  # Calculate proportions per cell-cycle stage
  cc_df <- plyr::ddply(seurat_object@meta.data, .(group), function(x) {
    out_df <- data.frame(table(x$Phase))
    return(out_df)
  })
  
  colnames(cc_df) <- c("group", "stage", "count")
  cc_df[is.na(cc_df)] <- 0
  
  cc_df <- cc_df %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(prop = count / sum(count)) %>%
    data.frame()
  
  # Generate bar plot
  p <- ggplot(data = cc_df, aes(x = group, y = prop, fill = stage)) +
    geom_col(colour = "black", alpha = 0.6, width = 0.95) +
    scale_fill_manual(values = colorRampPalette(
      brewer.pal("Spectral", n = 11)[c(3, 8, 10)])(3)) +
    theme_minimal() +
    xlab("\nGroup") +
    ylab("Proportion\n") + 
    theme(
      text = element_text(size = 22),
      panel.grid = element_blank(),
      axis.ticks.x = element_line()
    )
  
  # Adjust x-axis text for certain groupings
  if (group_diff == "sample") {
    p <- p + theme(
      axis.text.x = element_text(angle = 90, hjust = 1)
    )
  }
  
  # Save the plot
  ggsave2(file = paste0("CC_Stage_by_", group_diff, ".pdf"),
          plot = p,
          height = 8, 
          width = 7,
          bg = "white")
  
  # Return to the original directory
  setwd("../")
}


#' Average Expression Correlation
#'
#' This function calculates the correlation of average expression between groups in a Seurat object 
#' and visualises the correlation as a heatmap.
#'
#' @param seurat_object A Seurat object containing the data.
#' @param group_diff The grouping variable (e.g., "sample", "type") for which to calculate correlations.
#' @param output_dir Directory to save the output heatmap (default: "non_int_analysis_SCT").
#'
#' @return Saves a heatmap showing the correlation of average expression between groups as a PDF file.
#'
#' @examples
#' scRNA_Avg_Exp_Corr_fun(
#'   seurat_object = seurat_comb_SCT,
#'   group_diff = "sample"
#' )
scRNA_Avg_Exp_Corr_fun <- function(seurat_object, group_diff, output_dir = "non_int_analysis_SCT") {
  
  # Create or move to the output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  setwd(output_dir)
  
  # Calculate average expression
  avg_exp_df <- AverageExpression(
    object = seurat_object,
    assays = "SCT",
    group.by = group_diff
  )[[1]] %>% data.frame()
  
  # Compute correlations
  cor_exp_df <- as.data.frame(cor(avg_exp_df))
  
  # Clean row and column names
  rownames(cor_exp_df) <- sub("X", "", rownames(cor_exp_df))
  colnames(cor_exp_df) <- sub("X", "", colnames(cor_exp_df))
  
  # Define colour palette
  colvec <- rev(c("#9E0142", "#D53E4F", "#F46D43", "#E6F598", "#ABDDA4", "#66C2A5", 
                  "#3288BD", "#5E4FA2"))
  
  # Save the heatmap
  pdf(file = paste0("Avg_Exp_Heatmap_by_", group_diff, ".pdf"), 
      units = "cm", res = 300, height = 12, width = 14)
  
  p <- Heatmap(cor_exp_df, col = colvec, name = "Corr")
  draw(p)
  dev.off()
  
  # Return to the original directory
  setwd("../")
}


#' Differential Gene Expression Analysis
#'
#' This function performs differential gene expression (DE) analysis on a Seurat object, 
#' either for all pairwise group comparisons or for each group independently. The results 
#' are saved as CSV files and visualised using ridge plots and UMAP feature plots.
#'
#' @param seurat_object A Seurat object to be analysed.
#' @param group_diff The grouping variable (e.g., "sample", "type") for DE analysis.
#' @param all_pairwise Logical, whether to perform all pairwise comparisons (default: FALSE).
#' @param output_dir Directory to save DE results and plots (default: "non_int_analysis_SCT").
#'
#' @return Saves DE results, ridge plots, and UMAP plots to files and returns the updated Seurat object.
#'
#' @examples
#' scRNA_DE_fun(
#'   seurat_object = seurat_comb_SCT,
#'   group_diff = "sample",
#'   all_pairwise = TRUE
#' )
scRNA_DE_fun <- function(seurat_object, group_diff, all_pairwise = FALSE, output_dir = "non_int_analysis_SCT") {
  
  # Create or move to the main output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  setwd(output_dir)
  
  # Create or move to the DE Analysis directory
  if (!dir.exists("DE_Analysis")) {
    dir.create("DE_Analysis")
  }
  setwd("DE_Analysis")
  
  # Set group identity
  Idents(seurat_object) <- group_diff
  seurat_object@meta.data$group <- Idents(seurat_object)
  
  # Perform DE analysis
  if (all_pairwise) {
    pw_type_combs <- t(combn(unique(seurat_object@meta.data$group), 2))
    de_markers <- sapply(1:nrow(pw_type_combs), function(x) {
      FindMarkers(seurat_object, 
                  ident.1 = pw_type_combs[x, 1], 
                  ident.2 = pw_type_combs[x, 2], 
                  logfc.threshold = 0.5, 
                  min.pct = 0.5, 
                  min.diff.pct = 0.2)
    }, simplify = FALSE)
  } else {
    de_markers <- sapply(unique(seurat_object@meta.data$group), function(x) {
      FindMarkers(seurat_object, 
                  ident.1 = x, 
                  logfc.threshold = 0.5, 
                  min.pct = 0.5, 
                  min.diff.pct = 0.2)
    }, simplify = FALSE)
  }
  
  # Name the DE results
  if (all_pairwise) {
    names(de_markers) <- sapply(1:length(de_markers), function(x) {
      paste(pw_type_combs[x, ], collapse = "_")
    })
  } else {
    names(de_markers) <- unique(seurat_object@meta.data$group)
  }
  
  # Filter for non-empty results
  de_markers <- de_markers[sapply(de_markers, function(x) { nrow(x) >= 1 })]
  
  # Stop if no DE genes are found
  if (length(de_markers) == 0) {
    cat("\nNo DE genes found. Ending DE analysis for this object.\n")
    setwd("../../")
    return(seurat_object)
  }
  
  # Sort and extract top DE genes
  de_markers <- lapply(de_markers, function(x) {
    x <- x[order(abs(x$avg_log2FC), decreasing = TRUE), ]
    return(x)
  })
  
  de_genes <- sapply(de_markers, rownames, USE.NAMES = TRUE, simplify = FALSE)
  de_genes <- sapply(de_genes, function(x) {
    if (length(x) > 20) {
      return(x[1:20])
    } else {
      return(x)
    }
  }, USE.NAMES = TRUE, simplify = FALSE)
  
  # Generate ridge plots
  Ridge_plot_list <- lapply(seq_along(de_genes), function(x) {
    if (all_pairwise) {
      RidgePlot(seurat_object, 
                idents = strsplit(names(de_genes)[x], split = "_")[[1]], 
                features = de_genes[[x]], 
                ncol = 4)
    } else {
      RidgePlot(seurat_object, features = de_genes[[x]], ncol = 4)
    }
  })
  
  for (i in seq_along(Ridge_plot_list)) {
    p <- Ridge_plot_list[[i]]
    comp_name <- names(de_genes)[[i]]
    ggsave2(filename = paste0("DE_Exp_ridge_comp_", comp_name, ".pdf"),
            plot = p, 
            height = ceiling(length(de_genes[[i]]) / 4) * 3, 
            width = 12)
  }
  
  # Generate UMAP plots
  UMAP_by_DE_list <- lapply(seq_along(de_genes), function(x) {
    FeaturePlot(seurat_object, features = de_genes[[x]], split.by = "type", by.col = FALSE)
  })
  
  for (i in seq_along(UMAP_by_DE_list)) {
    p <- UMAP_by_DE_list[[i]]
    comp_name <- names(de_genes)[[i]]
    ggsave2(filename = paste0("DE_Exp_UMAPS_comp_", comp_name, ".pdf"),
            plot = p,
            height = 6, width = length(de_genes[[i]]) * 2)
  }
  
  # Combine DE markers into a single data frame
  for (i in seq_along(de_markers)) {
    de_markers[[i]]["gene"] <- rownames(de_markers[[i]])
    de_markers[[i]]["comp"] <- names(de_markers)[[i]]
  }
  
  de_genes_comp_df <- bind_rows(de_markers)
  write.csv(de_genes_comp_df, 
            file = paste0("DE_genes_df_", group_diff, ".csv"),
            row.names = FALSE)
  
  # Return to the original directory
  setwd("../../")
  
  return(seurat_object)
}


#' Gene Set Enrichment Analysis
#'
#' This function performs Gene Set Enrichment (GSE) analysis using `fgsea` and `msigdbr`. 
#' It supports both single-group and pairwise comparisons. Results are saved as CSV files 
#' and visualised with dot plots.
#'
#' @param seurat_object A Seurat object for GSE analysis.
#' @param group_diff The grouping variable (e.g., "sample", "type") for GSE analysis.
#' @param all_pairwise Logical, whether to perform all pairwise comparisons (default: FALSE).
#' @param output_dir Directory to save GSE results and plots (default: "non_int_analysis_SCT").
#'
#' @return Saves GSE results as CSV files and dot plots as PDFs in the specified output directory.
#'
#' @examples
#' scRNA_GSE_Analysis_fun(
#'   seurat_object = seurat_comb_SCT,
#'   group_diff = "sample",
#'   all_pairwise = TRUE
#' )
scRNA_GSE_Analysis_fun <- function(seurat_object, group_diff, all_pairwise = FALSE, output_dir = "non_int_analysis_SCT") {
  
  # Create or move to the main output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  setwd(output_dir)
  
  # Create or move to the GSE analysis directory
  if (!dir.exists("GSE_analysis")) {
    dir.create("GSE_analysis")
  }
  setwd("GSE_analysis")
  
  # Set group identity
  Idents(seurat_object) <- group_diff
  seurat_object@meta.data$group <- Idents(seurat_object)
  
  # Generate pairwise combinations if needed
  if (all_pairwise) {
    pw_type_combs <- t(combn(unique(seurat_object@meta.data$group), 2))
  }
  
  # Helper function for GSE analysis
  get_GSE <- function(seur_obj, grp_dff, pw) {
    GSE.genes <- wilcoxauc(seur_obj, grp_dff)
    m_df <- msigdbr(species = "Homo sapiens", category = "H")
    fgsea_sets <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
    
    # Process groups for GSE
    group.GSE.genes <- lapply(unique(GSE.genes$group), function(x) {
      out_df <- GSE.genes %>%
        dplyr::filter(group == x) %>%
        arrange(desc(auc)) %>%
        dplyr::select(feature, auc)
      tibble::deframe(out_df)
    })
    
    names(group.GSE.genes) <- unique(GSE.genes$group)
    
    # Run GSE with fgsea
    group.GSE.fgseaRes <- lapply(group.GSE.genes, function(x) {
      fgseaRes <- fgseaMultilevel(fgsea_sets, stats = x)
      fgseaRes$pathway <- sub("HALLMARK_", "", fgseaRes$pathway)
      fgseaRes
    })
    
    # Add group names to results
    for (i in seq_along(unique(GSE.genes$group))) {
      group.GSE.fgseaRes[[i]]$group <- unique(GSE.genes$group)[[i]]
    }
    
    # Combine results
    GSE.fgseaRes <- bind_rows(group.GSE.fgseaRes) %>% select(-leadingEdge)
    
    # Save results
    if (pw) {
      pw_comp <- paste0(unique(GSE.fgseaRes$group), collapse = "_")
      write.csv(GSE.fgseaRes, file = paste0("GSE_fgseaRes_pw_", pw_comp, ".csv"), row.names = FALSE)
    } else {
      write.csv(GSE.fgseaRes, file = paste0("GSE_fgseaRes_", group_diff, ".csv"), row.names = FALSE)
    }
    
    # Select top pathways for visualisation
    sub.GSE.fgseaRes <- lapply(group.GSE.fgseaRes, function(x) {
      x %>% arrange(desc(NES)) %>% head(n = 20)
    })
    
    bind_rows(sub.GSE.fgseaRes)
  }
  
  # Perform GSE analysis
  if (all_pairwise) {
    for (i in 1:nrow(pw_type_combs)) {
      seurat_object_sub <- subset(seurat_object, subset = (group == pw_type_combs[i, 1]) | (group == pw_type_combs[i, 2]))
      GSE_df <- get_GSE(seurat_object_sub, group_diff, TRUE)
      
      # Plot results
      p <- ggplot(subset(GSE_df, padj < 0.05), aes(x = group, y = reorder(pathway, NES), fill = NES, size = -log10(padj))) +
        geom_point(shape = 21, colour = "black") +
        scale_fill_viridis_c() +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        ylab("") +
        xlab("")
      
      pw_comp <- paste0(pw_type_combs[i, ], collapse = "_")
      ggsave2(filename = paste0("GSE_by_pw_", pw_comp, ".pdf"), plot = p, height = 8, width = 6 + length(unique(GSE_df$group)) * 0.2, bg = "white")
    }
  } else {
    GSE_df <- get_GSE(seurat_object, group_diff, FALSE)
    
    # Plot results
    p <- ggplot(subset(GSE_df, padj < 0.05), aes(x = group, y = reorder(pathway, NES), fill = NES, size = -log10(padj))) +
      geom_point(shape = 21, colour = "black") +
      scale_fill_viridis_c() +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ylab("") +
      xlab("")
    
    ggsave2(filename = paste0("GSE_by_", group_diff, ".pdf"), plot = p, height = 8, width = 6 + length(unique(GSE_df$group)) * 0.2, bg = "white")
  }
  
  # Return to the original directory
  setwd("../../")
}
