
################################################################################
# edgeR Differential Expression analysis R script. 
# Freddie Whiting
# Jan 2025
################################################################################

# This script performs differential expression analysis using the edgeR package.
# It includes functions to categorize samples into groups, filter comparison 
# dataframes based on specified criteria, and load necessary libraries for 
# the analysis. Specific comparisons are made by identifying genes that are 
# differentially expressed in one sample group comparison but are not 
# differentially expressed in another. For example, we isolate genes associated
# with resistance by finding genes that are differentially expressed in the 
# parental (POT) vs drug-treatment (DT) replicates, but that are not 
# differentially expressed in the parental vs acute drug exposure (DS)
# replicates. The script also includes functions to perform gene set enrichment
# analysis using fgsea and msigdbr.

################################################################################

rm(list = ls())

# Load necessary libraries:

library(plyr)
library(dplyr)
library(purrr)
library(ggplot2)
library(cowplot)
library(limma)
library(MAST)
library(tidyr)
library(msigdbr)
library(fgsea)
library(presto)
library(RColorBrewer)
library(knitr)
library(edgeR)
library(Matrix.utils)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(eulerr)

setwd("path/to/Barcode_Analysis/QR_10X_scRNA")

################################################################################

###########
# Functions
###########

#' @title Sample to Group Assignment
#' @description Assigns a group type to each sample based on naming patterns.
#' @param samp Character vector of sample names.
#' @return Character vector of group names.
#' @examples
#' samp_to_group(c("CO_1", "DS_2", "POT_3"))
samp_to_group <- function(samp){
  if(length(grep("CO", samp)) == 1){
    group <- "CO"
  } else if(length(grep("DS", samp)) == 1){
    group <- "DS"
  } else {
    group <- sub("(+{2,3})_.*", replacement = "\\1", samp)
  }
  return(group)
}
samp_to_group <- Vectorize(samp_to_group)


#' @title Filter Comparison Dataframe
#' @description Adds a filter column to a dataframe based on specified criteria.
#' @param comp_df Dataframe containing comparison data.
#' @param fdr_lim Numeric. Threshold for FDR.
#' @param pval_lim Numeric. Threshold for p-value.
#' @param lfc_lim Numeric. Threshold for absolute log fold-change.
#' @param fdr_lt Logical. If TRUE, filters for FDR < fdr_lim, otherwise FDR > fdr_lim.
#' @return A dataframe with an additional filter column.
#' @examples
#' filter_comp_fun(df, 0.05, 0.01, 1.0, TRUE)
filter_comp_fun <- function(comp_df, 
                            fdr_lim, pval_lim, lfc_lim, 
                            fdr_lt = TRUE){
  if(fdr_lt){
    comp_df <- comp_df %>%
      mutate(filter = case_when(
        FDR < fdr_lim & PValue < pval_lim & abs(logFC) > lfc_lim ~ 1,
        TRUE ~ 0
      ))
  } else {
    comp_df <- comp_df %>%
      mutate(filter = case_when(
        FDR > fdr_lim ~ 1,
        TRUE ~ 0
      ))
  }
  comp_df$filter <- as.factor(comp_df$filter)
  return(comp_df)
}


#' @title Aggregate Counts for Comparison
#' @description Aggregates counts from two groups in a comparison dataframe.
#' @param comp_df Dataframe with comparison data.
#' @param count_agg_df Dataframe with aggregated count data.
#' @return A dataframe with aggregated counts and their sum for the comparison.
#' @examples
#' comp_agg_fun(comparison_df, aggregated_counts_df)
comp_agg_fun <- function(comp_df, count_agg_df){
  
  comp_elements <- (comp_df$comparison %>% 
                      unique() %>% 
                      strsplit(split = " - "))[[1]] %>%
    as.list()
  
  # Get the sum of the RNA counts per gene for the two comparison elements:
  
  comp_agg_counts <- lapply(comp_elements, 
                            function(x){
                              count_agg_df %>%
                                subset(group == x) %>%
                                dplyr::select(comp_df$gene) %>%
                                apply(2, sum)
                            })
  
  for(i in seq_along(comp_agg_counts)){
    comp_agg_counts[[i]] <- data.frame(gene = names(comp_agg_counts[[i]]),
                                       count = comp_agg_counts[[i]],
                                       group = comp_elements[[i]])
    colnames(comp_agg_counts[[i]])[2] <- paste0("count_", i)
    colnames(comp_agg_counts[[i]])[3] <- paste0("group_", i)
    rownames(comp_agg_counts[[i]]) <- NULL
  }
  
  comp_agg_counts <- purrr::reduce(comp_agg_counts, full_join)
  comp_agg_counts$count_sum <- comp_agg_counts$count_1 + comp_agg_counts$count_2
  return(comp_agg_counts)
  
}


#' @title Filter and Combine Single Comparison
#' @description Filters and combines data for a single comparison, adding relevant columns.
#' @param comp Dataframe with comparison data.
#' @param count_agg_df Dataframe with aggregated counts.
#' @param pval_lim Numeric. P-value threshold.
#' @param lfc_lim Numeric. Log fold-change threshold.
#' @param fdr_lim Numeric. FDR threshold.
#' @param fdr_lt Logical. If TRUE, filters for FDR < fdr_lim, otherwise FDR > fdr_lim.
#' @return A dataframe with filtered and combined data.
#' @examples
#' fil_comb_fun(comp, count_agg_df, 0.05, 1.0, 0.05, TRUE)
fil_comb_fun <- function(comp, count_agg_df, 
                         pval_lim, lfc_lim,
                         fdr_lim, 
                         fdr_lt){
  comp_df <- filter_comp_fun(comp, fdr_lim, pval_lim, lfc_lim, fdr_lt=fdr_lt)
  comp_df_agg <- comp_agg_fun(comp_df, count_agg_df)

  comp_sub <- comp_df %>%
    dplyr::select(gene, logFC, PValue, filter) %>%
    rename_with(~ paste0("comp_1_", .x), -gene)

  comp_sub <- full_join(comp_sub, comp_df_agg, by = "gene")
  comp_sub$comp_zc_filter <- as.numeric(
    comp_sub$count_1 > 0 & comp_sub$count_2 > 0
  )
  return(comp_sub)
}


#' @title Filter and Combine Two Comparisons
#' @description Filters and combines data for two comparisons, adding relevant columns.
#' @param comp_1 Dataframe for the first comparison.
#' @param comp_2 Dataframe for the second comparison.
#' @param count_agg_df Dataframe with aggregated counts.
#' @param pval_lim Numeric. P-value threshold for filtering.
#' @param lfc_lim Numeric. Log fold-change threshold for filtering.
#' @param fdr_lim_1 Numeric. FDR threshold for the first comparison.
#' @param fdr_lim_2 Numeric. FDR threshold for the second comparison.
#' @param fdr_lt_1 Logical. If TRUE, filters for FDR < fdr_lim_1; otherwise FDR > fdr_lim_1.
#' @param fdr_lt_2 Logical. If TRUE, filters for FDR < fdr_lim_2; otherwise FDR > fdr_lim_2.
#' @return A dataframe with filtered and combined data for both comparisons, including additional filter and count-related columns.
#' @examples
#' fil_comb_comps_fun(comp_1, comp_2, count_agg_df, 0.05, 1.0, 0.05, 0.1, TRUE, FALSE)
fil_comb_comps_fun <- function(comp_1, comp_2, count_agg_df,
                               pval_lim, lfc_lim, 
                               fdr_lim_1, fdr_lim_2, 
                               fdr_lt_1, fdr_lt_2) {
  comp_1 <- filter_comp_fun(comp_1, fdr_lim_1, pval_lim, lfc_lim, fdr_lt = fdr_lt_1)
  comp_2 <- filter_comp_fun(comp_2, fdr_lim_2, pval_lim, lfc_lim, fdr_lt = fdr_lt_2)

  comp_1_agg <- comp_agg_fun(comp_1, count_agg_df)
  comp_2_agg <- comp_agg_fun(comp_2, count_agg_df)

  colnames(comp_2_agg) <- gsub("_1", "_3", colnames(comp_2_agg))
  colnames(comp_2_agg) <- gsub("_2", "_4", colnames(comp_2_agg))

  comp_1_sub <- comp_1 %>%
    dplyr::select(gene, logFC, PValue, filter) %>%
    rename_with(~ paste0("comp_1_", .x), -gene)

  comp_2_sub <- comp_2 %>%
    dplyr::select(gene, logFC, PValue, filter) %>%
    rename_with(~ paste0("comp_2_", .x), -gene)

  comp_sub <- full_join(comp_1_sub, comp_2_sub, by = "gene")
  comp_sub <- full_join(comp_sub, comp_1_agg, by = "gene")
  comp_sub <- full_join(comp_sub, comp_2_agg, by = "gene")

  comp_sub$comp_1_zc_filter <- as.numeric(
    comp_sub$count_1 > 0 & comp_sub$count_2 > 0
  )
  comp_sub$comp_2_zc_filter <- as.numeric(
    comp_sub$count_3 > 0 & comp_sub$count_4 > 0
  )
  comp_sub$comp_both_filter <- (comp_sub$comp_1_filter == "1") + 
                               (comp_sub$comp_2_filter == "1") == 2

  return(comp_sub)
}


# Pathway analysis comparison plot functions:
#############################################

#' @title Pathway Analysis Plot for FGSEA
#' @description Generates a comparative pathway enrichment plot using FGSEA results for two datasets.
#' @param fgsea_df1 Dataframe of FGSEA results for the first dataset.
#' @param fgsea_df2 Dataframe of FGSEA results for the second dataset.
#' @param samp_1 Label for the first dataset.
#' @param samp_2 Label for the second dataset.
#' @param pval_lim Numeric. P-value threshold for significance.
#' @param comp_name Character. Title for the comparison plot.
#' @return A ggplot object representing the pathway comparison.
#' @examples
#' comp_fgsea_plot(fgsea_df1, fgsea_df2, "Sample A", "Sample B", 0.05, "Comparison")
comp_fgsea_plot <- function(fgsea_df1, fgsea_df2,
                            samp_1, samp_2,
                            pval_lim,
                            comp_name) {
  fgsea_df1$pathway <- gsub("HALLMARK_", "", fgsea_df1$pathway)
  fgsea_df2$pathway <- gsub("HALLMARK_", "", fgsea_df2$pathway)

  fgsea_df1 <- fgsea_df1 %>% arrange(padj)
  fgsea_df2 <- fgsea_df2 %>% arrange(padj)
  pathway_order <- unique(c(fgsea_df1$pathway, fgsea_df2$pathway))
  fgsea_df1$pathway <- factor(fgsea_df1$pathway, levels = pathway_order)
  fgsea_df2$pathway <- factor(fgsea_df2$pathway, levels = pathway_order)

  if (nrow(subset(fgsea_df1, padj <= 0.05)) +
      nrow(subset(fgsea_df2, padj <= 0.05)) == 0) {
    p <- ggplot() + geom_blank()
  } else {
    p <- ggplot() +
      geom_point(data = subset(fgsea_df1, padj <= pval_lim),
                 aes(x = samp_1,
                     y = pathway,
                     fill = NES,
                     size = -log(padj)),
                 shape = 21, colour = "black", alpha = 0.8) +
      geom_point(data = subset(fgsea_df2, padj <= pval_lim),
                 aes(x = samp_2,
                     y = pathway,
                     fill = NES,
                     size = -log(padj)),
                 shape = 21, colour = "black", alpha = 0.8) +
      scale_fill_distiller(palette = "Spectral",
                           direction = -1) +
      scale_size_continuous(range = c(0.1, 10)) +
      xlab("") +
      ylab("") +
      ggtitle(comp_name) +
      theme_minimal()
  }

  return(p)
}


#' @title Pathway Analysis Plot for KEGG
#' @description Generates a comparative pathway enrichment plot using KEGG results for two datasets.
#' @param kegg_df1 Dataframe of KEGG pathway results for the first dataset.
#' @param kegg_df2 Dataframe of KEGG pathway results for the second dataset.
#' @param samp_1 Label for the first dataset.
#' @param samp_2 Label for the second dataset.
#' @param pval_lim Numeric. P-value threshold for filtering pathways.
#' @param comp_name Character. Title for the comparison plot.
#' @return A ggplot object representing the KEGG pathway comparison.
#' @examples
#' comp_kegg_plot(kegg_df1, kegg_df2, "Sample A", "Sample B", 0.05, "Comparison")
comp_kegg_plot <- function(kegg_df1, kegg_df2, 
                           samp_1, samp_2, 
                           pval_lim,
                           comp_name) {
  kegg_df1 <- topKEGG(kegg_df1, p.value = pval_lim, number = 20)
  kegg_df2 <- topKEGG(kegg_df2, p.value = pval_lim, number = 20)

  pathway_order <- full_join(kegg_df1 %>%
                               dplyr::select(Pathway, P.ENTREZID) %>%
                               dplyr::rename(p1 = P.ENTREZID),
                             kegg_df2 %>%
                               dplyr::select(Pathway, P.ENTREZID) %>%
                               dplyr::rename(p2 = P.ENTREZID)) %>%
    arrange(replace_na(p1, Inf),
            replace_na(p2, Inf)) %>%
    dplyr::pull(Pathway) %>% rev()

  kegg_df1$Pathway <- factor(kegg_df1$Pathway, levels = pathway_order)
  kegg_df2$Pathway <- factor(kegg_df2$Pathway, levels = pathway_order)

  if (nrow(kegg_df1) + nrow(kegg_df2) == 0) {
    p <- ggplot() + geom_blank()
  } else {
    p <- ggplot() + 
      geom_point(data = kegg_df1, 
                 aes(x = samp_1,
                     y = Pathway,
                     fill = ENTREZID,
                     size = -log(P.ENTREZID)),
                 shape = 21, colour = "black", alpha = 0.8) + 
      geom_point(data = kegg_df2, 
                 aes(x = samp_2,
                     y = Pathway,
                     fill = ENTREZID,
                     size = -log(P.ENTREZID)),
                 shape = 21, colour = "black", alpha = 0.8) + 
      scale_fill_viridis_c(name = "n genes") +
      scale_size_continuous(range = c(0.1, 10),
                            name = "-log(pval)") +
      xlab("") + 
      ylab("") + 
      ggtitle(comp_name) +
      theme_minimal()
  }
  
  return(p)
}


#' @title Pathway Analysis Plot for Gene Ontology (GO)
#' @description Generates a comparative pathway enrichment plot using GO results for two datasets.
#' @param go_df1 Dataframe of GO enrichment results for the first dataset.
#' @param go_df2 Dataframe of GO enrichment results for the second dataset.
#' @param samp_1 Label for the first dataset.
#' @param samp_2 Label for the second dataset.
#' @param pval_lim Numeric. P-value threshold for filtering pathways.
#' @param comp_name Character. Title for the comparison plot.
#' @return A ggplot object representing the GO term comparison.
#' @examples
#' comp_go_plot(go_df1, go_df2, "Sample A", "Sample B", 0.05, "Comparison")
comp_go_plot <- function(go_df1, go_df2, 
                         samp_1, samp_2, 
                         pval_lim,
                         comp_name) {
  go_df1 <- topGO(go_df1, p.value = pval_lim, number = 20)
  go_df2 <- topGO(go_df2, p.value = pval_lim, number = 20)

  pathway_order <- full_join(go_df1 %>%
                               dplyr::select(Term, P.ENTREZID) %>%
                               dplyr::rename(p1 = P.ENTREZID),
                             go_df2 %>%
                               dplyr::select(Term, P.ENTREZID) %>%
                               dplyr::rename(p2 = P.ENTREZID)) %>%
    arrange(replace_na(p1, Inf),
            replace_na(p2, Inf)) %>%
    dplyr::pull(Term) %>% rev()

  go_df1$Term <- factor(go_df1$Term, levels = pathway_order)
  go_df2$Term <- factor(go_df2$Term, levels = pathway_order)

  if (nrow(go_df1) + nrow(go_df2) == 0) {
    p <- ggplot() + geom_blank()
  } else {
    p <- ggplot() + 
      geom_point(data = go_df1, 
                 aes(x = samp_1,
                     y = Term,
                     fill = ENTREZID,
                     size = -log(P.ENTREZID)),
                 shape = 21, colour = "black", alpha = 0.8) + 
      geom_point(data = go_df2, 
                 aes(x = samp_2,
                     y = Term,
                     fill = ENTREZID,
                     size = -log(P.ENTREZID)),
                 shape = 21, colour = "black", alpha = 0.8) + 
      scale_fill_viridis_c(name = "n genes") +
      scale_size_continuous(range = c(0.1, 10),
                            name = "-log(pval)") +
      xlab("") + 
      ylab("") + 
      ggtitle(comp_name) +
      theme_minimal()
  }

  return(p)
}


#' @title Differential Expression Analysis with edgeR
#' @description Performs differential expression (DE) analysis using the edgeR package on sample count data.
#' @param samp_pb Matrix of sample counts, where rows are genes and columns are samples.
#' @param run_regexp Character. A regular expression to determine which contrasts to perform based on the dataset.
#' @return A dataframe containing DE results for all contrasts, including logFC, p-values, FDR, and the associated comparisons.
#' @examples
#' edgeR_DE_analysis(sample_counts_matrix, "HCTbc")
edgeR_DE_analysis <- function(samp_pb, run_regexp) {
  
  # Convert the count matrix to a DGEList object
  dge_samp <- DGEList(t(samp_pb))
  
  # Assign groups based on sample names
  samples <- rownames(dge_samp$samples)
  groups <- samp_to_group(samples)
  dge_samp$samples$group <- as.factor(groups)
  
  # Design matrix for the DE analysis
  design <- model.matrix(~0 + group, data = dge_samp$samples)
  colnames(design) <- levels(dge_samp$samples$group)
  
  # Filter out lowly expressed genes
  keep <- filterByExpr(dge_samp, 
                       min.count = 10,
                       min.total.count = 40)
  dge_samp <- dge_samp[keep, , keep.lib.sizes = FALSE]
  
  # Estimate dispersion
  dge_samp <- estimateDisp(dge_samp, design)
  
  # Plot biological coefficient of variation (BCV)
  plotBCV(dge_samp)
  
  # Fit the generalized linear model with quasi-likelihood
  fit <- glmQLFit(dge_samp, design, robust = TRUE)
  
  # Plot quasi-likelihood dispersions
  plotQLDisp(fit)
  
  # Define contrasts based on the dataset identifier
  if (grepl("HCTbc", run_regexp)) {
    my.contrasts <- makeContrasts("DT3 - POTs",
                                  "DT4 - POTs",
                                  "DS - POTs",
                                  "DT3 - DS",
                                  "DT4 - DS",
                                  "DT3 - DT4",
                                  "DT4 - DT3",
                                  "CO - POTs",
                                  levels = design)
  } else if (grepl("SW6bc", run_regexp)) {
    my.contrasts <- makeContrasts("DT1 - POTs",
                                  "DT3 - POTs",
                                  "DS - POTs",
                                  "DT1 - DS",
                                  "DT3 - DS",
                                  "DT1 - DT3",
                                  "DT3 - DT1",
                                  "CO - POTs",
                                  levels = design)
  }
  
  # Perform the quasi-likelihood F-test for each contrast
  qlf <- list()
  for (i in seq_along(colnames(my.contrasts))) {
    qlf[[i]] <- glmQLFTest(fit, contrast = my.contrasts[, i])
    qlf[[i]]$comparison <- colnames(my.contrasts)[i]
  }
  
  # Compile DE results into a single dataframe
  all_top_tags_df <- lapply(qlf, function(x) {
    top_df <- topTags(x, n = Inf)$table
    top_df$comparison <- x$comparison
    top_df$gene <- rownames(top_df)
    rownames(top_df) <- NULL
    return(top_df)
  }) %>% bind_rows()
  
  return(all_top_tags_df)
}


#################
# Hallmark Subset
#################

hm_subset_df <- read.csv("path/to/Barcode_Analysis/QR_10X_scRNA/subset_gsea_hallmark_pathways_Jan_2024.csv")

##########
# Plotting
##########

# Custom colour palettes for plotting: 

Hco_pal <- colorRampPalette(brewer.pal("YlGnBu", n = 9)[3:6])(4)
Hdt_pal <- colorRampPalette(brewer.pal("RdPu", n = 9)[6:9])(4)
Sco_pal <- colorRampPalette(brewer.pal("PuBu", n = 9)[3:6])(4)
Sdt_pal <- colorRampPalette(brewer.pal("OrRd", n = 9)[6:9])(4)

Hds_pal <- colorRampPalette(brewer.pal("YlOrBr", n = 9)[3:5])(3)
Sds_pal <- colorRampPalette(brewer.pal("YlOrRd", n = 9)[3:5])(3)

HPOT_col <- colorRampPalette(brewer.pal("Greens", n = 9)[4])(1)
SPOT_col <- colorRampPalette(brewer.pal("Greens", n = 9)[7])(1)

################################################################################

##################################
# Differential Expression Analysis
##################################

# Read in the corresponding aggregated count information for each cell line:

QR_HCTbc_samp_pb <- read.csv("QR_HCTbc/seurat_analysis/QR_HCTbc_edgeR_DE_sample_agg_counts.csv",
                             row.names=1)
QR_SW6bc_samp_pb <- read.csv("QR_SW6bc/seurat_analysis/QR_SW6bc_edgeR_DE_sample_agg_counts.csv",
                             row.names=1)

# Run edgeR:

HCT_DE_df <- edgeR_DE_analysis(QR_HCTbc_samp_pb, "QR_HCTbc")
SW6_DE_df <- edgeR_DE_analysis(QR_SW6bc_samp_pb, "QR_SW6bc")

################################################################################

# Make a note of the entrez gene IDs:

QR_HCTbc_entrez_genes <- select(org.Hs.eg.db,
                                keys = colnames(QR_HCTbc_samp_pb),
                                columns = "ENTREZID",
                                keytype = "SYMBOL") %>%
  subset(!is.na(ENTREZID))

QR_SW6bc_entrez_genes <- select(org.Hs.eg.db,
                                keys = colnames(QR_SW6bc_samp_pb),
                                columns = "ENTREZID",
                                keytype = "SYMBOL") %>%
  subset(!is.na(ENTREZID))


##########
# QR_SW6bc
##########

setwd("path/to/Barcode_Analysis/QR_10X_scRNA")
setwd("QR_SW6bc/")
setwd("seurat_analysis/")

SW6_samp_agg_df <- read.csv("QR_SW6bc_edgeR_DE_sample_agg_counts.csv") %>%
  dplyr::rename(samp = X)

# Check gene names are the same in each: 

setdiff(colnames(SW6_samp_agg_df),
        unique(SW6_DE_df$gene))

SW6_samp_agg_df$group <- samp_to_group(SW6_samp_agg_df$samp)

# Split by the comparisons: 

SW6_DE_df_split <- SW6_DE_df %>%
  group_split(comparison)

# Choose two comparisons to explore simultaneously: 

lapply(SW6_DE_df_split, function(x){unique(x$comparison)})

comparison_order <- lapply(SW6_DE_df_split, function(x){unique(x$comparison)})


# Go through and collect the dataframes that make all the specific comparisons
# between samples we're interested in: 

# DT1 Resistance Genes
######################

# Simultaneously checks if genes are DE in DT1 vs POT and not DE in DS v POT:

S_DT1_res_df <- fil_comb_comps_fun(SW6_DE_df_split[[5]],
                                   SW6_DE_df_split[[2]],
                                   SW6_samp_agg_df,
                                   0.10, 1.0,
                                   0.05, 0.1,
                                   T, F)

# We want to include a score that we can rank the genes by that captures the 
# two comparisons (different in DT vs POT whilst little change in DS vs POT)
# whilst retaining all genes so that we can perform GSEA. 

S_DT1_res_df$rank_score <- S_DT1_res_df$comp_1_logFC * (1 / (abs(S_DT1_res_df$comp_2_logFC)+1))

S_DT1_res_df <- S_DT1_res_df %>%
  arrange(rank_score)


# DT1 Acute Exposure Genes
##########################

# Simultaneously checks if genes are DE in DT3 vs POT and DE in DS v POT:

S_DT1_acu_df <- fil_comb_comps_fun(SW6_DE_df_split[[5]],
                                   SW6_DE_df_split[[2]],
                                   SW6_samp_agg_df,
                                   0.10, 2.0, 
                                   0.05, 0.05, 
                                   T, T)

# We now want a rank that simply incorporates both logFCs

S_DT1_acu_df$rank_score <- S_DT1_acu_df$comp_1_logFC * abs(S_DT1_acu_df$comp_2_logFC+1)

S_DT1_acu_df <- S_DT1_acu_df %>%
  arrange(rank_score)
       

# DS (vs DT1) Persister Genes
#############################

# Simultaneously checks if genes are DE in DS vs POT and not DE in DT1 v POT:

S_DT1_per_df <- fil_comb_comps_fun(SW6_DE_df_split[[2]],
                                   SW6_DE_df_split[[5]],
                                   SW6_samp_agg_df,
                                   0.10, 1.0,
                                   0.05, 0.1,
                                   T, F)

# We want to include a score that we can rank the genes by that captures the 
# two comparisons (different in DT vs POT whilst little change in DS vs POT)
# whilst retaining all genes so that we can perform GSEA. 

S_DT1_per_df$rank_score <- S_DT1_per_df$comp_1_logFC * (1 / (abs(S_DT1_per_df$comp_2_logFC)+1))

S_DT1_per_df <- S_DT1_per_df %>%
  arrange(rank_score)


# DT3 Resistance Genes
######################

# Simultaneously checks if genes are DE in DT3 vs POT and not DE in DS v POT:

S_DT3_res_df <- fil_comb_comps_fun(SW6_DE_df_split[[8]],
                                   SW6_DE_df_split[[2]],
                                   SW6_samp_agg_df,
                                   0.10, 1.0,
                                   0.05, 0.1,
                                   T, F)

# We want to include a score that we can rank the genes by that captures the 
# two comparisons (different in DT vs POT whilst little change in DS vs POT)
# whilst retaining all genes so that we can perform GSEA. 

S_DT3_res_df$rank_score <- S_DT3_res_df$comp_1_logFC * (1 / (abs(S_DT3_res_df$comp_2_logFC)+1))

S_DT3_res_df <- S_DT3_res_df %>%
  arrange(rank_score)


# DT3 Acute Exposure Genes
##########################

# Simultaneously checks if genes are DE in DT3 vs POT and DE in DS v POT:

S_DT3_acu_df <- fil_comb_comps_fun(SW6_DE_df_split[[8]],
                                   SW6_DE_df_split[[2]],
                                   SW6_samp_agg_df,
                                   0.10, 2.0, 
                                   0.05, 0.05, 
                                   T, T)

# We now want a rank that simply incorporates both logFCs

S_DT3_acu_df$rank_score <- S_DT3_acu_df$comp_1_logFC * abs(S_DT3_acu_df$comp_2_logFC+1)

S_DT3_acu_df <- S_DT3_acu_df %>%
  arrange(rank_score)


# DS (vs DT3) Persister Genes
#############################

# Simultaneously checks if genes are DE in DS vs POT and not DE in DT1 v POT:

S_DT3_per_df <- fil_comb_comps_fun(SW6_DE_df_split[[2]],
                                   SW6_DE_df_split[[8]],
                                   SW6_samp_agg_df,
                                   0.10, 1.0,
                                   0.05, 0.1,
                                   T, F)

# We want to include a score that we can rank the genes by that captures the 
# two comparisons (different in DT vs POT whilst little change in DS vs POT)
# whilst retaining all genes so that we can perform GSEA. 

S_DT3_per_df$rank_score <- S_DT3_per_df$comp_1_logFC * (1 / (abs(S_DT3_per_df$comp_2_logFC)+1))

S_DT3_per_df <- S_DT3_per_df %>%
  arrange(rank_score)

# Number of genes per comparison:
#################################

S_DT1_nres_up <- subset(S_DT1_res_df, comp_1_logFC > 0 & comp_both_filter==T) %>% 
  nrow()
S_DT1_nres_dn <- subset(S_DT1_res_df, comp_1_logFC < 0 & comp_both_filter==T) %>% 
  nrow()

S_DT1_nacu_up <- subset(S_DT1_acu_df, comp_1_logFC > 0 & comp_both_filter==T) %>% 
  nrow()
S_DT1_nacu_dn <- subset(S_DT1_acu_df, comp_1_logFC < 0 & comp_both_filter==T) %>% 
  nrow()

S_DT1_nper_up <- subset(S_DT1_per_df, comp_1_logFC > 0 & comp_both_filter==T) %>% 
  nrow()
S_DT1_nper_dn <- subset(S_DT1_per_df, comp_1_logFC < 0 & comp_both_filter==T) %>% 
  nrow()

S_DT3_nres_up <- subset(S_DT3_res_df, comp_1_logFC > 0 & comp_both_filter==T) %>% 
  nrow()
S_DT3_nres_dn <- subset(S_DT3_res_df, comp_1_logFC < 0 & comp_both_filter==T) %>% 
  nrow()

S_DT3_nacu_up <- subset(S_DT3_acu_df, comp_1_logFC > 0 & comp_both_filter==T) %>% 
  nrow()
S_DT3_nacu_dn <- subset(S_DT3_acu_df, comp_1_logFC < 0 & comp_both_filter==T) %>% 
  nrow()

S_DT3_nper_up <- subset(S_DT3_per_df, comp_1_logFC > 0 & comp_both_filter==T) %>% 
  nrow()
S_DT3_nper_dn <- subset(S_DT3_per_df, comp_1_logFC < 0 & comp_both_filter==T) %>% 
  nrow()


S_DT_comp_ngenes_df <- data.frame(samp = c(rep("S_DT1", 6), rep("S_DT3", 6)),
                                  comp = rep(c("res", "res", "acu", "acu", "per", "per"), 2),
                                  dir = rep(c("up", "down"), 6),
                                  n_genes = c(S_DT1_nres_up, -S_DT1_nres_dn,
                                              S_DT1_nacu_up, -S_DT1_nacu_dn,
                                              S_DT1_nper_up, -S_DT1_nper_dn,
                                              S_DT3_nres_up, -S_DT3_nres_dn,
                                              S_DT3_nacu_up, -S_DT3_nacu_dn,
                                              S_DT3_nper_up, -S_DT3_nper_dn))


# Gene Set Enrichment
#####################

# Get the msigdbr hallmark gene sets: 

m_df<- msigdbr(species = "Homo sapiens", category = "H")

fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

# Subset by given hallmark pathways:

fgsea_sets <- fgsea_sets[names(fgsea_sets) %in% hm_subset_df$pathway]

# Rank the genes:
#################

# DT1 Res Ranks:
S_DT1_res_ranks <- S_DT1_res_df %>%
  subset(gene %in% QR_SW6bc_entrez_genes$SYMBOL) %>% 
  dplyr::arrange(rank_score) %>%
  dplyr::select(gene, rank_score) %>%
  tibble::deframe()

# DT1 Acu Ranks: 
S_DT1_acu_ranks <- S_DT1_acu_df %>%
  subset(gene %in% QR_SW6bc_entrez_genes$SYMBOL) %>% 
  dplyr::arrange(rank_score) %>%
  dplyr::select(gene, rank_score) %>%
  tibble::deframe()

# DT1 Per Ranks:
S_DT1_per_ranks <- S_DT1_per_df %>%
  subset(gene %in% QR_SW6bc_entrez_genes$SYMBOL) %>% 
  dplyr::arrange(rank_score) %>%
  dplyr::select(gene, rank_score) %>%
  tibble::deframe()


# DT3 Res Ranks: 
S_DT3_res_ranks <- S_DT3_res_df %>%
  subset(gene %in% QR_SW6bc_entrez_genes$SYMBOL) %>% 
  dplyr::arrange(rank_score) %>%
  dplyr::select(gene, rank_score) %>%
  tibble::deframe()

# DT3 Acu Ranks: 
S_DT3_acu_ranks <- S_DT3_acu_df %>%
  subset(gene %in% QR_SW6bc_entrez_genes$SYMBOL) %>% 
  dplyr::arrange(desc(rank_score)) %>%
  dplyr::select(gene, rank_score) %>%
  tibble::deframe()

# DT3 Per Ranks: 
S_DT3_per_ranks <- S_DT3_per_df %>%
  subset(gene %in% QR_SW6bc_entrez_genes$SYMBOL) %>% 
  dplyr::arrange(desc(rank_score)) %>%
  dplyr::select(gene, rank_score) %>%
  tibble::deframe()

# FGSEA:
########

S_DT1_res_fgsea <- fgseaMultilevel(fgsea_sets, stats = S_DT1_res_ranks)
S_DT3_res_fgsea <- fgseaMultilevel(fgsea_sets, stats = S_DT3_res_ranks)

S_DT1_acu_fgsea <- fgseaMultilevel(fgsea_sets, stats = S_DT1_acu_ranks)
S_DT3_acu_fgsea <- fgseaMultilevel(fgsea_sets, stats = S_DT3_acu_ranks)

S_DT1_per_fgsea <- fgseaMultilevel(fgsea_sets, stats = S_DT1_per_ranks)
S_DT3_per_fgsea <- fgseaMultilevel(fgsea_sets, stats = S_DT3_per_ranks)


# KEGG
######

# DT1 Res KEGG terms:
S_DT1_res_kegg <- QR_SW6bc_entrez_genes %>%
  subset(SYMBOL %in% subset(S_DT1_res_df, comp_both_filter==T)$gene) %>%
  dplyr::select("ENTREZID") %>%
  as.vector() %>%
  kegga(species = "Hs")


# DT1 Acu KEGG terms:
S_DT1_acu_kegg <- QR_SW6bc_entrez_genes %>%
  subset(SYMBOL %in% subset(S_DT1_acu_df, comp_both_filter==T)$gene) %>%
  dplyr::select("ENTREZID") %>%
  as.vector() %>%
  kegga(species = "Hs")


# DT1 Per KEGG terms:
S_DT1_per_kegg <- QR_SW6bc_entrez_genes %>%
  subset(SYMBOL %in% subset(S_DT1_per_df, comp_both_filter==T)$gene) %>%
  dplyr::select("ENTREZID") %>%
  as.vector() %>%
  kegga(species = "Hs")


# DT3 Res KEGG terms:
S_DT3_res_kegg <- QR_SW6bc_entrez_genes %>%
  subset(SYMBOL %in% subset(S_DT3_res_df, comp_both_filter==T)$gene) %>%
  dplyr::select("ENTREZID") %>%
  as.vector() %>%
  kegga(species = "Hs")


# DT3 Acu KEGG terms:
S_DT3_acu_kegg <- QR_SW6bc_entrez_genes %>%
  subset(SYMBOL %in% subset(S_DT3_acu_df, comp_both_filter==T)$gene) %>%
  dplyr::select("ENTREZID") %>%
  as.vector() %>%
  kegga(species = "Hs")

topKEGG(S_DT3_acu_kegg)

# DT3 Per KEGG terms:
S_DT3_per_kegg <- QR_SW6bc_entrez_genes %>%
  subset(SYMBOL %in% subset(S_DT3_per_df, comp_both_filter==T)$gene) %>%
  dplyr::select("ENTREZID") %>%
  as.vector() %>%
  kegga(species = "Hs")



# GO
####

# DT1 Res GO terms:
S_DT1_res_goana <- QR_SW6bc_entrez_genes %>%
  subset(SYMBOL %in% subset(S_DT1_res_df, comp_both_filter==T)$gene) %>%
  dplyr::select("ENTREZID") %>%
  as.vector() %>%
  goana(species = "Hs")


# DT1 Acu GO terms:
S_DT1_acu_goana <- QR_SW6bc_entrez_genes %>%
  subset(SYMBOL %in% subset(S_DT1_acu_df, comp_both_filter==T)$gene) %>%
  dplyr::select("ENTREZID") %>%
  as.vector() %>%
  goana(species = "Hs")


# DT1 Per GO terms:
S_DT1_per_goana <- QR_SW6bc_entrez_genes %>%
  subset(SYMBOL %in% subset(S_DT1_per_df, comp_both_filter==T)$gene) %>%
  dplyr::select("ENTREZID") %>%
  as.vector() %>%
  goana(species = "Hs")


# DT3 Res GO terms:
S_DT3_res_goana <- QR_SW6bc_entrez_genes %>%
  subset(SYMBOL %in% subset(S_DT3_res_df, comp_both_filter==T)$gene) %>%
  dplyr::select("ENTREZID") %>%
  as.vector() %>%
  goana(species = "Hs")


# DT3 Acu GO terms:
S_DT3_acu_goana <- QR_SW6bc_entrez_genes %>%
  subset(SYMBOL %in% subset(S_DT3_acu_df, comp_both_filter==T)$gene) %>%
  dplyr::select("ENTREZID") %>%
  as.vector() %>%
  goana(species = "Hs")


# DT1 Per GO terms:
S_DT3_per_goana <- QR_SW6bc_entrez_genes %>%
  subset(SYMBOL %in% subset(S_DT3_per_df, comp_both_filter==T)$gene) %>%
  dplyr::select("ENTREZID") %>%
  as.vector() %>%
  goana(species = "Hs")



# DT1 vs DT3
############

# Get a list of genes that are 'resistant' genes across both replicates:

S_DT1_res_genes <- subset(S_DT1_res_df, comp_both_filter==T)$gene
S_DT3_res_genes <- subset(S_DT3_res_df, comp_both_filter==T)$gene

S_DT_res_genes <- unique(S_DT1_res_genes, S_DT3_res_genes)

# DT1 vs DT3
# Differentially expressed resistance associated genes:

S_DT_comp_dif_df <- fil_comb_fun(SW6_DE_df_split[[4]],
                                 SW6_samp_agg_df,
                                 0.10, 2.0, 
                                 0.05, T)

S_DT_dif_genes <- S_DT_comp_dif_df %>%
  subset(gene %in% S_DT_res_genes) %>%
  subset(comp_1_filter==1) %>%
  dplyr::select(gene) %>% as.vector()


S_DT_comp_dif_df <- S_DT_comp_dif_df %>%
  mutate(dif_exp = gene %in% S_DT_dif_genes$gene)



# DT1 vs DT3
# Shared expressed resistance associated genes:

S_DT_comp_shr_df <- fil_comb_fun(SW6_DE_df_split[[4]],
                                 SW6_samp_agg_df,
                                 0.10, 1.0, 
                                 0.10, F)

S_DT_shr_genes <- S_DT_comp_shr_df %>%
  subset(gene %in% S_DT_res_genes) %>%
  subset(comp_1_filter==1) %>%
  dplyr::select(gene) %>% as.vector()

length(S_DT_shr_genes$gene)

S_DT_comp_shr_df <- S_DT_comp_shr_df %>%
  mutate(shr_exp = gene %in% S_DT_shr_genes$gene)


# Combine DT1 and DT3 rank for GSEA:

S_DT1_res_rank_df <- S_DT1_res_df %>% dplyr::select(gene, rank_score) %>% dplyr::rename(DT1_rank_score=rank_score)

S_DT3_res_rank_df <- S_DT3_res_df %>% dplyr::select(gene, rank_score) %>% dplyr::rename(DT3_rank_score=rank_score)

S_DT_res_rank_df <- full_join(S_DT1_res_rank_df, 
                              S_DT3_res_rank_df) %>%
  dplyr::mutate(rank_sum = DT1_rank_score+DT3_rank_score,
                rank_diff = DT1_rank_score-DT3_rank_score)


# Rank by difference in ranks: 
S_DT_dif_ranks <- S_DT_res_rank_df %>%
  subset(gene %in% QR_SW6bc_entrez_genes$SYMBOL) %>% 
  dplyr::arrange(desc(rank_diff)) %>%
  dplyr::select(gene, rank_diff) %>%
  tibble::deframe()

# Rank by sum of ranks:
S_DT_sum_ranks <- S_DT_res_rank_df %>%
  subset(gene %in% QR_SW6bc_entrez_genes$SYMBOL) %>% 
  dplyr::arrange(desc(rank_sum)) %>%
  dplyr::select(gene, rank_sum) %>%
  tibble::deframe()


# FGSEA:
########

S_DT_dif_fgsea <- fgseaMultilevel(fgsea_sets, stats = S_DT_dif_ranks)
S_DT_sum_fgsea <- fgseaMultilevel(fgsea_sets, stats = S_DT_sum_ranks)


# KEGG
######

# DT Dif Res KEGG terms:
S_DT_dif_res_kegg <- QR_SW6bc_entrez_genes %>%
  subset(SYMBOL %in% S_DT_dif_genes$gene) %>%
  dplyr::select("ENTREZID") %>%
  as.vector() %>%
  kegga(species = "Hs")

# DT Shr Acu KEGG terms:
S_DT_shr_res_kegg <- QR_SW6bc_entrez_genes %>%
  subset(SYMBOL %in% S_DT_shr_genes$gene) %>%
  dplyr::select("ENTREZID") %>%
  as.vector() %>%
  kegga(species = "Hs")

# GO
####

# DT Diff res GO terms:
S_DT_dif_goana <- QR_SW6bc_entrez_genes %>%
  subset(SYMBOL %in% S_DT_dif_genes$gene) %>%
  dplyr::select("ENTREZID") %>%
  as.vector() %>%
  goana(species = "Hs")

# DT Shr res GO terms:
S_DT_shr_goana <- QR_SW6bc_entrez_genes %>%
  subset(SYMBOL %in% S_DT_shr_genes$gene) %>%
  dplyr::select("ENTREZID") %>%
  as.vector() %>%
  goana(species = "Hs")

##########
# QR_HCTbc
##########

setwd("path/to/Barcode_Analysis/QR_10X_scRNA")
setwd("QR_HCTbc/")
setwd("seurat_analysis/")

HCT_samp_agg_df <- read.csv("QR_HCTbc_edgeR_DE_sample_agg_counts.csv") %>%
  dplyr::rename(samp = X)

# Check gene names are the same in each: 

HCT_samp_agg_df$group <- samp_to_group(HCT_samp_agg_df$samp)

# Split by the comparisons: 

HCT_DE_df_split <- HCT_DE_df %>%
  group_split(comparison)

# Choose two comparisons to explore simultaneously: 

comparison_order <- lapply(HCT_DE_df_split, function(x){unique(x$comparison)})


# DT3 Resistance Genes
######################

# Simultaneously checks if genes are DE in DT3 vs POT and not DE in DS v POT:

H_DT3_res_df <- fil_comb_comps_fun(HCT_DE_df_split[[5]],
                                   HCT_DE_df_split[[2]],
                                   HCT_samp_agg_df,
                                   0.10, 1.0,
                                   0.05, 0.1,
                                   T, F)

# We want to include a score that we can rank the genes by that captures the 
# two comparisons (different in DT vs POT whilst little change in DS vs POT)
# whilst retaining all genes so that we can perform GSEA. 

H_DT3_res_df$rank_score <- H_DT3_res_df$comp_1_logFC * (1 / (abs(H_DT3_res_df$comp_2_logFC)+1))

H_DT3_res_df <- H_DT3_res_df %>%
  arrange(rank_score)


# DT3 Acute Exposure Genes
##########################

# Simultaneously checks if genes are DE in DT3 vs POT and DE in DS v POT:

H_DT3_acu_df <- fil_comb_comps_fun(HCT_DE_df_split[[5]],
                                   HCT_DE_df_split[[2]],
                                   HCT_samp_agg_df,
                                   0.10, 2.0, 
                                   0.05, 0.05, 
                                   T, T)

# We now want a rank that simply incorporates both logFCs

H_DT3_acu_df$rank_score <- H_DT3_acu_df$comp_1_logFC * abs(H_DT3_acu_df$comp_2_logFC+1)

H_DT3_acu_df <- H_DT3_acu_df %>%
  arrange(rank_score)
       

# DS (vs DT3) Persister Genes
#############################

# Simultaneously checks if genes are DE in DS vs POT and not DE in DT3 v POT:

H_DT3_per_df <- fil_comb_comps_fun(HCT_DE_df_split[[2]],
                                   HCT_DE_df_split[[5]],
                                   HCT_samp_agg_df,
                                   0.10, 1.0,
                                   0.05, 0.1,
                                   T, F)

# We want to include a score that we can rank the genes by that captures the 
# two comparisons (different in DT vs POT whilst little change in DS vs POT)
# whilst retaining all genes so that we can perform GSEA. 

H_DT3_per_df$rank_score <- H_DT3_per_df$comp_1_logFC * (1 / (abs(H_DT3_per_df$comp_2_logFC)+1))

H_DT3_per_df <- H_DT3_per_df %>%
  arrange(rank_score)


# DT4 Resistance Genes
######################

# Simultaneously checks if genes are DE in DT4 vs POT and not DE in DS v POT:

H_DT4_res_df <- fil_comb_comps_fun(HCT_DE_df_split[[8]],
                                   HCT_DE_df_split[[2]],
                                   HCT_samp_agg_df,
                                   0.10, 1.0,
                                   0.05, 0.1,
                                   T, F)

# We want to include a score that we can rank the genes by that captures the 
# two comparisons (different in DT vs POT whilst little change in DS vs POT)
# whilst retaining all genes so that we can perform GSEA. 

H_DT4_res_df$rank_score <- H_DT4_res_df$comp_1_logFC * (1 / (abs(H_DT4_res_df$comp_2_logFC)+1))

H_DT4_res_df <- H_DT4_res_df %>%
  arrange(rank_score)


# DT4 Acute Exposure Genes
##########################

# Simultaneously checks if genes are DE in DT4 vs POT and DE in DS v POT:

H_DT4_acu_df <- fil_comb_comps_fun(HCT_DE_df_split[[8]],
                                   HCT_DE_df_split[[2]],
                                   HCT_samp_agg_df,
                                   0.10, 2.0, 
                                   0.05, 0.05, 
                                   T, T)

# We now want a rank that simply incorporates both logFCs

H_DT4_acu_df$rank_score <- H_DT4_acu_df$comp_1_logFC * abs(H_DT4_acu_df$comp_2_logFC+1)

H_DT4_acu_df <- H_DT4_acu_df %>%
  arrange(rank_score)


# DS (vs DT4) Persister Genes
#############################

# Simultaneously checks if genes are DE in DS vs POT and not DE in DT3 v POT:

H_DT4_per_df <- fil_comb_comps_fun(HCT_DE_df_split[[2]],
                                   HCT_DE_df_split[[8]],
                                   HCT_samp_agg_df,
                                   0.10, 1.0,
                                   0.05, 0.1,
                                   T, F)

# We want to include a score that we can rank the genes by that captures the 
# two comparisons (different in DT vs POT whilst little change in DS vs POT)
# whilst retaining all genes so that we can perform GSEA. 

H_DT4_per_df$rank_score <- H_DT4_per_df$comp_1_logFC * (1 / (abs(H_DT4_per_df$comp_2_logFC)+1))

H_DT4_per_df <- H_DT4_per_df %>%
  arrange(rank_score)


# Number of genes per comparison:
#################################

H_DT3_nres_up <- subset(H_DT3_res_df, comp_1_logFC > 0 & comp_both_filter==T) %>% 
  nrow()
H_DT3_nres_dn <- subset(H_DT3_res_df, comp_1_logFC < 0 & comp_both_filter==T) %>% 
  nrow()

H_DT3_nacu_up <- subset(H_DT3_acu_df, comp_1_logFC > 0 & comp_both_filter==T) %>% 
  nrow()
H_DT3_nacu_dn <- subset(H_DT3_acu_df, comp_1_logFC < 0 & comp_both_filter==T) %>% 
  nrow()

H_DT3_nper_up <- subset(H_DT3_per_df, comp_1_logFC > 0 & comp_both_filter==T) %>% 
  nrow()
H_DT3_nper_dn <- subset(H_DT3_per_df, comp_1_logFC < 0 & comp_both_filter==T) %>% 
  nrow()

H_DT4_nres_up <- subset(H_DT4_res_df, comp_1_logFC > 0 & comp_both_filter==T) %>% 
  nrow()
H_DT4_nres_dn <- subset(H_DT4_res_df, comp_1_logFC < 0 & comp_both_filter==T) %>% 
  nrow()

H_DT4_nacu_up <- subset(H_DT4_acu_df, comp_1_logFC > 0 & comp_both_filter==T) %>% 
  nrow()
H_DT4_nacu_dn <- subset(H_DT4_acu_df, comp_1_logFC < 0 & comp_both_filter==T) %>% 
  nrow()

H_DT4_nper_up <- subset(H_DT4_per_df, comp_1_logFC > 0 & comp_both_filter==T) %>% 
  nrow()
H_DT4_nper_dn <- subset(H_DT4_per_df, comp_1_logFC < 0 & comp_both_filter==T) %>% 
  nrow()


H_DT_comp_ngenes_df <- data.frame(samp = c(rep("H_DT3", 6), rep("H_DT4", 6)),
                                  comp = rep(c("res", "res", "acu", "acu", "per", "per"), 2),
                                  dir = rep(c("up", "down"), 6),
                                  n_genes = c(H_DT3_nres_up, -H_DT3_nres_dn,
                                              H_DT3_nacu_up, -H_DT3_nacu_dn,
                                              H_DT3_nper_up, -H_DT3_nper_dn,
                                              H_DT4_nres_up, -H_DT4_nres_dn,
                                              H_DT4_nacu_up, -H_DT4_nacu_dn,
                                              H_DT4_nper_up, -H_DT4_nper_dn))


# Gene Set Enrichment
#####################

# Get the msigdbr hallmark gene sets: 

m_df<- msigdbr(species = "Homo sapiens", category = "H")

fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

# Subset by given hallmark pathways:

fgsea_sets <- fgsea_sets[names(fgsea_sets) %in% hm_subset_df$pathway]

# Rank the genes:
#################

# DT1 Res Ranks:
H_DT3_res_ranks <- H_DT3_res_df %>%
  subset(gene %in% QR_HCTbc_entrez_genes$SYMBOL) %>% 
  dplyr::arrange(rank_score) %>%
  dplyr::select(gene, rank_score) %>%
  tibble::deframe()

# DT1 Acu Ranks: 
H_DT3_acu_ranks <- H_DT3_acu_df %>%
  subset(gene %in% QR_HCTbc_entrez_genes$SYMBOL) %>% 
  dplyr::arrange(rank_score) %>%
  dplyr::select(gene, rank_score) %>%
  tibble::deframe()

# DT1 Per Ranks:
H_DT3_per_ranks <- H_DT3_per_df %>%
  subset(gene %in% QR_HCTbc_entrez_genes$SYMBOL) %>% 
  dplyr::arrange(rank_score) %>%
  dplyr::select(gene, rank_score) %>%
  tibble::deframe()


# DT3 Res Ranks: 
H_DT4_res_ranks <- H_DT4_res_df %>%
  subset(gene %in% QR_HCTbc_entrez_genes$SYMBOL) %>% 
  dplyr::arrange(rank_score) %>%
  dplyr::select(gene, rank_score) %>%
  tibble::deframe()

# DT3 Acu Ranks: 
H_DT4_acu_ranks <- H_DT4_acu_df %>%
  subset(gene %in% QR_HCTbc_entrez_genes$SYMBOL) %>% 
  dplyr::arrange(desc(rank_score)) %>%
  dplyr::select(gene, rank_score) %>%
  tibble::deframe()

# DT3 Per Ranks: 
H_DT4_per_ranks <- H_DT4_per_df %>%
  subset(gene %in% QR_HCTbc_entrez_genes$SYMBOL) %>% 
  dplyr::arrange(desc(rank_score)) %>%
  dplyr::select(gene, rank_score) %>%
  tibble::deframe()

# FGSEA:
########

H_DT3_res_fgsea <- fgseaMultilevel(fgsea_sets, stats = H_DT3_res_ranks)
H_DT4_res_fgsea <- fgseaMultilevel(fgsea_sets, stats = H_DT4_res_ranks)

H_DT3_acu_fgsea <- fgseaMultilevel(fgsea_sets, stats = H_DT3_acu_ranks)
H_DT4_acu_fgsea <- fgseaMultilevel(fgsea_sets, stats = H_DT4_acu_ranks)

H_DT3_per_fgsea <- fgseaMultilevel(fgsea_sets, stats = H_DT3_per_ranks)
H_DT4_per_fgsea <- fgseaMultilevel(fgsea_sets, stats = H_DT4_per_ranks)


# KEGG
######

# DT1 Res KEGG terms:
H_DT3_res_kegg <- QR_HCTbc_entrez_genes %>%
  subset(SYMBOL %in% subset(H_DT3_res_df, comp_both_filter==T)$gene) %>%
  dplyr::select("ENTREZID") %>%
  as.vector() %>%
  kegga(species = "Hs")

# DT1 Acu KEGG terms:
H_DT3_acu_kegg <- QR_HCTbc_entrez_genes %>%
  subset(SYMBOL %in% subset(H_DT3_acu_df, comp_both_filter==T)$gene) %>%
  dplyr::select("ENTREZID") %>%
  as.vector() %>%
  kegga(species = "Hs")

# DT1 Per KEGG terms:
H_DT3_per_kegg <- QR_HCTbc_entrez_genes %>%
  subset(SYMBOL %in% subset(H_DT3_per_df, comp_both_filter==T)$gene) %>%
  dplyr::select("ENTREZID") %>%
  as.vector() %>%
  kegga(species = "Hs")

# DT3 Res KEGG terms:
H_DT4_res_kegg <- QR_HCTbc_entrez_genes %>%
  subset(SYMBOL %in% subset(H_DT4_res_df, comp_both_filter==T)$gene) %>%
  dplyr::select("ENTREZID") %>%
  as.vector() %>%
  kegga(species = "Hs")

# DT3 Acu KEGG terms:
H_DT4_acu_kegg <- QR_HCTbc_entrez_genes %>%
  subset(SYMBOL %in% subset(H_DT4_acu_df, comp_both_filter==T)$gene) %>%
  dplyr::select("ENTREZID") %>%
  as.vector() %>%
  kegga(species = "Hs")

# DT3 Per KEGG terms:
H_DT4_per_kegg <- QR_HCTbc_entrez_genes %>%
  subset(SYMBOL %in% subset(H_DT4_per_df, comp_both_filter==T)$gene) %>%
  dplyr::select("ENTREZID") %>%
  as.vector() %>%
  kegga(species = "Hs")


# GO
####

# DT3 Res GO terms:
H_DT3_res_goana <- QR_HCTbc_entrez_genes %>%
  subset(SYMBOL %in% subset(H_DT3_res_df, comp_both_filter==T)$gene) %>%
  dplyr::select("ENTREZID") %>%
  as.vector() %>%
  goana(species = "Hs")

# DT3 Acu GO terms:
H_DT3_acu_goana <- QR_HCTbc_entrez_genes %>%
  subset(SYMBOL %in% subset(H_DT3_acu_df, comp_both_filter==T)$gene) %>%
  dplyr::select("ENTREZID") %>%
  as.vector() %>%
  goana(species = "Hs")

# DT3 Per GO terms:
H_DT3_per_goana <- QR_HCTbc_entrez_genes %>%
  subset(SYMBOL %in% subset(H_DT3_per_df, comp_both_filter==T)$gene) %>%
  dplyr::select("ENTREZID") %>%
  as.vector() %>%
  goana(species = "Hs")

# DT4 Res GO terms:
H_DT4_res_goana <- QR_HCTbc_entrez_genes %>%
  subset(SYMBOL %in% subset(H_DT4_res_df, comp_both_filter==T)$gene) %>%
  dplyr::select("ENTREZID") %>%
  as.vector() %>%
  goana(species = "Hs")

# DT4 Acu GO terms:
H_DT4_acu_goana <- QR_HCTbc_entrez_genes %>%
  subset(SYMBOL %in% subset(H_DT4_acu_df, comp_both_filter==T)$gene) %>%
  dplyr::select("ENTREZID") %>%
  as.vector() %>%
  goana(species = "Hs")

# DT4 Per GO terms:
H_DT4_per_goana <- QR_HCTbc_entrez_genes %>%
  subset(SYMBOL %in% subset(H_DT4_per_df, comp_both_filter==T)$gene) %>%
  dplyr::select("ENTREZID") %>%
  as.vector() %>%
  goana(species = "Hs")


# DT3 vs DT4
############

# Get a list of genes that are 'resistant' genes across both replicates:

H_DT3_res_genes <- subset(H_DT3_res_df, comp_both_filter==T)$gene
H_DT4_res_genes <- subset(H_DT4_res_df, comp_both_filter==T)$gene

H_DT_res_genes <- unique(H_DT3_res_genes, H_DT4_res_genes)

# DT3 vs DT4
# Differentially expressed resistance associated genes:

H_DT_comp_dif_df <- fil_comb_fun(HCT_DE_df_split[[4]],
                                 HCT_samp_agg_df,
                                 0.10, 2.0, 
                                 0.05, T)

H_DT_dif_genes <- H_DT_comp_dif_df %>%
  subset(gene %in% H_DT_res_genes) %>%
  subset(comp_1_filter==1) %>%
  dplyr::select(gene) %>% as.vector()

H_DT_comp_dif_df <- H_DT_comp_dif_df %>%
  mutate(dif_exp = gene %in% H_DT_dif_genes$gene)


# DT3 vs DT4
# Shared expressed resistance associated genes:

H_DT_comp_shr_df <- fil_comb_fun(HCT_DE_df_split[[4]],
                                 HCT_samp_agg_df,
                                 0.10, 1.0, 
                                 0.10, F)

H_DT_shr_genes <- H_DT_comp_shr_df %>%
  subset(gene %in% H_DT_res_genes) %>%
  subset(comp_1_filter==1) %>%
  dplyr::select(gene) %>% as.vector()


H_DT_comp_shr_df <- H_DT_comp_shr_df %>%
  mutate(shr_exp = gene %in% H_DT_shr_genes$gene)


# Combine DT3 and DT4 rank for GSEA:

H_DT3_res_rank_df <- H_DT3_res_df %>% dplyr::select(gene, rank_score) %>% dplyr::rename(DT1_rank_score=rank_score)

H_DT4_res_rank_df <- H_DT4_res_df %>% dplyr::select(gene, rank_score) %>% dplyr::rename(DT3_rank_score=rank_score)

H_DT_res_rank_df <- full_join(H_DT3_res_rank_df, 
                              H_DT4_res_rank_df) %>%
  dplyr::mutate(rank_sum = DT1_rank_score+DT3_rank_score,
                rank_diff = DT1_rank_score-DT3_rank_score)


# Rank by difference in ranks: 
H_DT_dif_ranks <- H_DT_res_rank_df %>%
  subset(gene %in% QR_HCTbc_entrez_genes$SYMBOL) %>% 
  dplyr::arrange(desc(rank_diff)) %>%
  dplyr::select(gene, rank_diff) %>%
  tibble::deframe()

# Rank by sum of ranks:
H_DT_sum_ranks <- H_DT_res_rank_df %>%
  subset(gene %in% QR_HCTbc_entrez_genes$SYMBOL) %>% 
  dplyr::arrange(desc(rank_sum)) %>%
  dplyr::select(gene, rank_sum) %>%
  tibble::deframe()


# FGSEA:
########

H_DT_dif_fgsea <- fgseaMultilevel(fgsea_sets, stats = H_DT_dif_ranks)
H_DT_sum_fgsea <- fgseaMultilevel(fgsea_sets, stats = H_DT_sum_ranks)


# KEGG
######

# DT Dif Res KEGG terms:
H_DT_dif_res_kegg <- QR_HCTbc_entrez_genes %>%
  subset(SYMBOL %in% H_DT_dif_genes$gene) %>%
  dplyr::select("ENTREZID") %>%
  as.vector() %>%
  kegga(species = "Hs")

# DT Shr Acu KEGG terms:
H_DT_shr_res_kegg <- QR_HCTbc_entrez_genes %>%
  subset(SYMBOL %in% H_DT_shr_genes$gene) %>%
  dplyr::select("ENTREZID") %>%
  as.vector() %>%
  kegga(species = "Hs")


# GO
####

# DT Diff res GO terms:
H_DT_dif_goana <- QR_HCTbc_entrez_genes %>%
  subset(SYMBOL %in% H_DT_dif_genes$gene) %>%
  dplyr::select("ENTREZID") %>%
  as.vector() %>%
  goana(species = "Hs")

# DT Shr res GO terms:
H_DT_shr_goana <- QR_HCTbc_entrez_genes %>%
  subset(SYMBOL %in% H_DT_shr_genes$gene) %>%
  dplyr::select("ENTREZID") %>%
  as.vector() %>%
  goana(species = "Hs")


################################################################################

#######################
# Save Comparison plots
#######################

#######
# SW6bc
#######

# GSEA
######

S_res_comp_fgsea_plot <- comp_fgsea_plot(S_DT1_res_fgsea,
                                         S_DT3_res_fgsea,
                                         "QR_SW6bc_DT1",
                                         "QR_SW6bc_DT3",
                                         0.05, "GSEA SW6bc: Resistance Genes")

S_acu_comp_fgsea_plot <- comp_fgsea_plot(S_DT1_acu_fgsea,
                                         S_DT3_acu_fgsea,
                                         "QR_SW6bc_DT1",
                                         "QR_SW6bc_DT3",
                                         0.05, "GSEA SW6bc: Acute Genes")

S_per_comp_fgsea_plot <- comp_fgsea_plot(S_DT1_per_fgsea,
                                         S_DT3_per_fgsea,
                                         "QR_SW6bc_DT1",
                                         "QR_SW6bc_DT3",
                                         0.05, "GSEA SW6bc: Persister Genes")

S_res_dif_shr_fgsea_plot <- comp_fgsea_plot(S_DT_dif_fgsea,
                                            S_DT_sum_fgsea,
                                            "QR_SW6bc Res Diff.",
                                            "QR_SW6bc Res Shared",
                                            0.05, "GSEA SW6bc: Resistant Diff. v Shared")

# KEGG
######

S_res_comp_kegg_plot <- comp_kegg_plot(S_DT1_res_kegg,
                                       S_DT3_res_kegg,
                                       "QR_SW6bc_DT1",
                                       "QR_SW6bc_DT3",
                                       0.05, "KEGG SW6bc: Resistance Genes")

S_acu_comp_kegg_plot <- comp_kegg_plot(S_DT1_acu_kegg,
                                       S_DT3_acu_kegg,
                                       "QR_SW6bc_DT1",
                                       "QR_SW6bc_DT3",
                                       0.05, "KEGG SW6bc: Acute Genes")

S_per_comp_kegg_plot <- comp_kegg_plot(S_DT1_per_kegg,
                                       S_DT3_per_kegg,
                                       "QR_SW6bc_DT1",
                                       "QR_SW6bc_DT3",
                                       0.05, "KEGG SW6bc: Persister Genes")

S_dif_shr_kegg_plot <- comp_kegg_plot(S_DT_dif_res_kegg,
                                      S_DT_shr_res_kegg,
                                      "QR_SW6bc Res Diff.",
                                      "QR_SW6bc Res Shared",
                                      0.05, "KEGG SW6bc: Resistant Diff. v Shared")

# GO
####

S_res_comp_go_plot <- comp_go_plot(S_DT1_res_goana,
                                   S_DT3_res_goana,
                                   "QR_SW6bc_DT1",
                                   "QR_SW6bc_DT3",
                                   0.05, "GO SW6bc: Resistance Genes")

S_acu_comp_go_plot <- comp_go_plot(S_DT1_acu_goana,
                                   S_DT3_acu_goana,
                                   "QR_SW6bc_DT1",
                                   "QR_SW6bc_DT3",
                                   0.05, "GO SW6bc: Acute Genes")

S_per_comp_go_plot <- comp_go_plot(S_DT1_per_goana,
                                   S_DT3_per_goana,
                                   "QR_SW6bc_DT1",
                                   "QR_SW6bc_DT3",
                                   0.05, "GO SW6bc: Persister Genes")

S_dif_shr_go_plot <- comp_go_plot(S_DT_dif_goana,
                                  S_DT_shr_goana,
                                  "QR_SW6bc Res Diff.",
                                  "QR_SW6bc Res Shared",
                                  0.05, "KEGG SW6bc: Resistant Diff. v Shared")

#######
# HCTbc
#######

# GSEA
######

H_res_comp_fgsea_plot <- comp_fgsea_plot(H_DT3_res_fgsea,
                                         H_DT4_res_fgsea,
                                         "QR_HCTbc_DT3",
                                         "QR_HCTbc_DT4",
                                         0.05, "GSEA HCTbc: Resistance Genes")

H_acu_comp_fgsea_plot <- comp_fgsea_plot(H_DT3_acu_fgsea,
                                         H_DT4_acu_fgsea,
                                         "QR_HCTbc_DT3",
                                         "QR_HCTbc_DT4",
                                         0.05, "GSEA HCTbc: Acute Genes")

H_per_comp_fgsea_plot <- comp_fgsea_plot(H_DT3_per_fgsea,
                                         H_DT4_per_fgsea,
                                         "QR_HCTbc_DT3",
                                         "QR_HCTbc_DT4",
                                         0.05, "GSEA HCTbc: Persister Genes")

H_res_dif_shr_fgsea_plot <- comp_fgsea_plot(H_DT_dif_fgsea,
                                            H_DT_sum_fgsea,
                                            "QR_HCTbc Res Diff.",
                                            "QR_HCTbc Res Shared",
                                            0.05, "GSEA HCTbc: Resistant Diff. v Shared")

# KEGG
######

H_res_comp_kegg_plot <- comp_kegg_plot(H_DT3_res_kegg,
                                       H_DT4_res_kegg,
                                       "QR_HCTbc_DT3",
                                       "QR_HCTbc_DT4",
                                       0.05, "KEGG HCTbc: Resistance Genes")

H_acu_comp_kegg_plot <- comp_kegg_plot(H_DT3_acu_kegg,
                                       H_DT4_acu_kegg,
                                       "QR_HCTbc_DT3",
                                       "QR_HCTbc_DT4",
                                       0.05, "KEGG HCTbc: Acute Genes")

H_per_comp_kegg_plot <- comp_kegg_plot(H_DT3_per_kegg,
                                       H_DT4_per_kegg,
                                       "QR_HCTbc_DT3",
                                       "QR_HCTbc_DT4",
                                       0.05, "KEGG HCTbc: Persister Genes")

H_dif_shr_kegg_plot <- comp_kegg_plot(H_DT_dif_res_kegg,
                                      H_DT_shr_res_kegg,
                                      "QR_HCTbc Res Diff.",
                                      "QR_HCTbc Res Shared",
                                      0.05, "KEGG HCTbc: Resistant Diff. v Shared")

# GO
####

H_res_comp_go_plot <- comp_go_plot(H_DT3_res_goana,
                                   H_DT4_res_goana,
                                   "QR_HCTbc_DT3",
                                   "QR_HCTbc_DT4",
                                   0.05, "GO HCTbc: Resistance Genes")

H_acu_comp_go_plot <- comp_go_plot(H_DT3_acu_goana,
                                   H_DT4_acu_goana,
                                   "QR_HCTbc_DT3",
                                   "QR_HCTbc_DT4",
                                   0.05, "GO HCTbc: Acute Genes")

H_per_comp_go_plot <- comp_go_plot(H_DT3_per_goana,
                                   H_DT4_per_goana,
                                   "QR_HCTbc_DT3",
                                   "QR_HCTbc_DT4",
                                   0.05, "GO HCTbc: Persister Genes")

H_dif_shr_go_plot <- comp_go_plot(H_DT_dif_goana,
                                  H_DT_shr_goana,
                                  "QR_HCTbc Res Diff.",
                                  "QR_HCTbc Res Shared",
                                  0.05, "KEGG HCTbc: Resistant Diff. v Shared")


plot_dir <- "path/to/Barcode_Analysis/QR_10X_scRNA/DE_Analysis_Plots"

if(dir.exists(plot_dir)==F){
  dir.create(plot_dir)
  setwd(plot_dir)
}else{
  setwd(plot_dir)
}

p_height <- 9
p_width <- 8

# PDFs
######

plots <- list(
    S_res_comp_fgsea_plot = S_res_comp_fgsea_plot,
    S_acu_comp_fgsea_plot = S_acu_comp_fgsea_plot,
    S_per_comp_fgsea_plot = S_per_comp_fgsea_plot,
    S_res_dif_shr_fgsea_plot = S_res_dif_shr_fgsea_plot,
    S_res_comp_kegg_plot = S_res_comp_kegg_plot,
    S_acu_comp_kegg_plot = S_acu_comp_kegg_plot,
    S_per_comp_kegg_plot = S_per_comp_kegg_plot,
    S_dif_shr_kegg_plot = S_dif_shr_kegg_plot,
    S_res_comp_go_plot = S_res_comp_go_plot,
    S_acu_comp_go_plot = S_acu_comp_go_plot,
    S_per_comp_go_plot = S_per_comp_go_plot,
    S_dif_shr_go_plot = S_dif_shr_go_plot,
    H_res_comp_fgsea_plot = H_res_comp_fgsea_plot,
    H_acu_comp_fgsea_plot = H_acu_comp_fgsea_plot,
    H_per_comp_fgsea_plot = H_per_comp_fgsea_plot,
    H_res_dif_shr_fgsea_plot = H_res_dif_shr_fgsea_plot,
    H_res_comp_kegg_plot = H_res_comp_kegg_plot,
    H_acu_comp_kegg_plot = H_acu_comp_kegg_plot,
    H_per_comp_kegg_plot = H_per_comp_kegg_plot,
    H_dif_shr_kegg_plot = H_dif_shr_kegg_plot,
    H_res_comp_go_plot = H_res_comp_go_plot,
    H_acu_comp_go_plot = H_acu_comp_go_plot,
    H_per_comp_go_plot = H_per_comp_go_plot,
    H_dif_shr_go_plot = H_dif_shr_go_plot
)

for (plot_name in names(plots)) {
    # Save as PDFs
    ggsave2(paste0(plot_name, ".pdf"),
                    plot = plots[[plot_name]],
                    height = p_height,
                    width = p_width,
                    device = cairo_pdf)
    # Save as JPGs
    ggsave(paste0(plot_name, ".jpg"),
                  plot = plots[[plot_name]],
                  height = p_height,
                  width = p_width)
                    
}

################################################################################

# Save Dataframes
#################

#######
# HCTbc
#######

plot_dir <- "path/to/Barcode_Analysis/QR_10X_scRNA/DE_Analysis_Plots"

if(dir.exists(plot_dir)==F){
  dir.create(plot_dir)
  setwd(plot_dir)
}else{
  setwd(plot_dir)
}

# DT resistance dataframes

write.csv(H_DT3_res_df,
          file = "QR_HCTbc_DT3_res_df.csv",
          row.names = F)

write.csv(H_DT4_res_df, 
          file = "QR_HCTbc_DT4_res_df.csv",
          row.names = F)

# DT acute exposure dataframes

write.csv(H_DT3_acu_df,
          file = "QR_HCTbc_DT3_acu_df.csv",
          row.names = F)

write.csv(H_DT4_acu_df,
          file = "QR_HCTbc_DT4_acu_df.csv",
          row.names = F)

# DT 'persistence' dataframes

write.csv(H_DT3_per_df,
          file = "QR_HCTbc_DT3_per_df.csv",
          row.names = F)

write.csv(H_DT4_per_df,
          file = "QR_HCTbc_DT4_per_df.csv",
          row.names = F)

# shr/dif resistance dataframes

write.csv(H_DT_comp_dif_df,
          file = "QR_HCTbc_DT_comp_dif_df.csv",
          row.names = F)

write.csv(H_DT_comp_shr_df, 
          file = "QR_HCTbc_DT_comp_shr_df.csv",
          row.names = F)


S_DT1_res_genes_sub <- S_DT1_res_df %>%
  subset(gene %in% QR_SW6bc_entrez_genes$SYMBOL) %>%
  subset(comp_both_filter==T) %>%
  pull(gene)

S_DT3_res_genes_sub <- S_DT3_res_df %>%
  subset(gene %in% QR_SW6bc_entrez_genes$SYMBOL) %>% 
  subset(comp_both_filter==T) %>%
  pull(gene)

H_DT3_res_genes_sub <- H_DT3_res_df %>%
  subset(gene %in% QR_HCTbc_entrez_genes$SYMBOL) %>%
  subset(comp_both_filter==T) %>%
  pull(gene)

H_DT4_res_genes_sub <- H_DT4_res_df %>%
  subset(gene %in% QR_HCTbc_entrez_genes$SYMBOL) %>% 
  subset(comp_both_filter==T) %>%
  pull(gene)



# SW6bc DT FGSEA Comparison
###########################

S_DT1_res_fgsea_df <- S_DT1_res_fgsea %>% 
  rename_with(~ if_else(.x != "pathway", paste0("DT1_", .x), .x))

S_DT3_res_fgsea_df <- S_DT3_res_fgsea %>% 
  rename_with(~ if_else(.x != "pathway", paste0("DT3_", .x), .x))

S_DT_res_fgsea_df <- dplyr::full_join(S_DT1_res_fgsea_df, S_DT3_res_fgsea_df) %>%
  subset(DT1_padj <= 0.05 | DT3_padj <= 0.05) %>% 
  dplyr::mutate(pathway = gsub("HALLMARK_", "", pathway))

S_DT_fgsea_comp_plot_1 <- ggplot() + 
  geom_point(data = S_DT_res_fgsea_df, aes(x = "QR_SW6bc_DT1",
                                           y = pathway,
                                           fill = DT1_NES,
                                           size = -log(DT1_padj)),
             shape = 21, colour = "black", alpha = 0.8) + 
  geom_point(data = S_DT_res_fgsea_df, aes(x = "QR_SW6bc_DT3",
                                         y = pathway,
                                         fill = DT3_NES,
                                         size = -log(DT3_padj)),
             shape = 21, colour = "black", alpha = 0.8) + 
  scale_fill_distiller(palette = "Spectral",
                       direction = -1,
                       name = "NES") + 
  scale_size_continuous(range = c(1.0, 12.0),
                        limits = c(0.0, 70.0),
                        name = "-log(p-value, adjusted)") +
  xlab("") + 
  ylab("") + 
  ggtitle("GSEA SW6bc: Resistance Genes") +
  theme_minimal()

S_DT_fgsea_comp_plot_2 <- ggplot(data = S_DT_res_fgsea_df, 
                                 aes(x = DT1_NES, y = DT3_NES,
                                     fill = DT1_NES)) +
  geom_point(shape = 21, colour = "black", size = 7, alpha = 0.6) +
  geom_text(aes(label = gsub("_", "\n", S_DT_res_fgsea_df$pathway)),
            size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_distiller(palette = "Spectral",
                     direction = -1,
                     name = "DT1 NES") + 
  xlab("\nDT1: NES") + 
  ylab("DT3: NES\n") + 
  theme_minimal() +
  theme(text = element_text(size = 7))


S_DT_fgsea_comp_plot <- cowplot::plot_grid(S_DT_fgsea_comp_plot_1,
                                           S_DT_fgsea_comp_plot_2,
                                           axis = c("lb"),
                                           align = "h",
                                           rel_widths = c(0.4, 0.5))

plot_dir <- "path/to/Barcode_Analysis/QR_10X_scRNA/DE_Analysis_Plots"

if(dir.exists(plot_dir)==F){
  dir.create(plot_dir)
  setwd(plot_dir)
}else{
  setwd(plot_dir)
}

p_height <- 8
p_width <- 15

ggsave2("QR_SW6bc_res_fgsea_cor_comp_plot.pdf",
        plot = S_DT_fgsea_comp_plot,
        height = p_height,
        width = p_width,
        device = cairo_pdf)

ggsave("QR_SW6bc_res_fgsea_cor_comp_plot.jpg",
        plot = S_DT_fgsea_comp_plot,
        height = p_height,
        width = p_width)


# SW6bc DS FGSEA Comparison
###########################

S_DT1_per_fgsea_df <- S_DT1_per_fgsea %>% 
  rename_with(~ if_else(.x != "pathway", paste0("DT1_", .x), .x))

S_DT3_per_fgsea_df <- S_DT3_per_fgsea %>% 
  rename_with(~ if_else(.x != "pathway", paste0("DT3_", .x), .x))

S_DT_per_fgsea_df <- dplyr::full_join(S_DT1_per_fgsea_df, S_DT3_per_fgsea_df) %>%
  subset(DT1_padj <= 0.05 | DT3_padj <= 0.05) %>% 
  dplyr::mutate(pathway = gsub("HALLMARK_", "", pathway))

S_DS_fgsea_comp_plot_1 <- ggplot() + 
  geom_point(data = S_DT_per_fgsea_df, aes(x = "QR_SW6bc_DS\n(vs DT1)",
                                           y = pathway,
                                           fill = DT1_NES,
                                           size = -log(DT1_padj)),
             shape = 21, colour = "black", alpha = 0.8) + 
  geom_point(data = S_DT_per_fgsea_df, aes(x = "QR_SW6bc_DS\n(vs DT3)",
                                         y = pathway,
                                         fill = DT3_NES,
                                         size = -log(DT3_padj)),
             shape = 21, colour = "black", alpha = 0.8) + 
  scale_fill_distiller(palette = "Spectral",
                       direction = -1,
                       name = "NES") + 
  scale_size_continuous(range = c(1.0, 12.0),
                        limits = c(0.0, 70.0),
                        name = "-log(p-value, adjusted)") +
  xlab("") + 
  ylab("") + 
  ggtitle("GSEA SW6bc: Persister Genes") +
  theme_minimal()

S_DS_fgsea_comp_plot_2 <- ggplot(data = S_DT_per_fgsea_df, 
                                 aes(x = DT1_NES, y = DT3_NES,
                                     fill = DT1_NES)) +
  geom_point(shape = 21, colour = "black", size = 7, alpha = 0.6) +
  geom_text(aes(label = gsub("_", "\n", S_DT_per_fgsea_df$pathway)),
            size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_distiller(palette = "Spectral",
                     direction = -1,
                     name = "DS (vs DT1)\nNES") + 
  xlab("\nDS (vs DT1): NES") + 
  ylab("DT3 (vs DT3): NES\n") + 
  theme_minimal()


S_DS_fgsea_comp_plot <- cowplot::plot_grid(S_DS_fgsea_comp_plot_1,
                                           S_DS_fgsea_comp_plot_2,
                                           axis = c("lb"),
                                           align = "h",
                                           rel_widths = c(0.5, 0.5))

plot_dir <- "path/to/Barcode_Analysis/QR_10X_scRNA/DE_Analysis_Plots"

if(dir.exists(plot_dir)==F){
  dir.create(plot_dir)
  setwd(plot_dir)
}else{
  setwd(plot_dir)
}

p_height <- 8
p_width <- 15

ggsave2("QR_SW6bc_per_fgsea_cor_comp_plot.pdf",
        plot = S_DS_fgsea_comp_plot,
        height = p_height,
        width = p_width,
        device = cairo_pdf)

ggsave2("QR_SW6bc_per_fgsea_cor_comp_plot.jpg",
        plot = S_DS_fgsea_comp_plot,
        height = p_height,
        width = p_width)


# HCTbc DT FGSEA Comparison
###########################

H_DT3_res_fgsea_df <- H_DT3_res_fgsea %>% 
  rename_with(~ if_else(.x != "pathway", paste0("DT3_", .x), .x))

H_DT4_res_fgsea_df <- H_DT4_res_fgsea %>% 
  rename_with(~ if_else(.x != "pathway", paste0("DT4_", .x), .x))

H_DT_res_fgsea_df <- dplyr::full_join(H_DT3_res_fgsea_df, H_DT4_res_fgsea_df) %>%
  subset(DT3_padj <= 0.05 | DT4_padj <= 0.05) %>% 
  dplyr::mutate(pathway = gsub("HALLMARK_", "", pathway))

H_DT_fgsea_comp_plot_1 <- ggplot() + 
  geom_point(data = H_DT_res_fgsea_df, aes(x = "QR_HCTbc_DT3",
                                           y = pathway,
                                           fill = DT3_NES,
                                           size = -log(DT3_padj)),
             shape = 21, colour = "black", alpha = 0.8) + 
  geom_point(data = H_DT_res_fgsea_df, aes(x = "QR_HCTbc_DT4",
                                         y = pathway,
                                         fill = DT4_NES,
                                         size = -log(DT4_padj)),
             shape = 21, colour = "black", alpha = 0.8) + 
  scale_fill_distiller(palette = "Spectral",
                       direction = -1,
                       name = "NES") + 
  scale_size_continuous(range = c(1.0, 12.0),
                        limits = c(0.0, 70.0),
                        name = "-log(p-value, adjusted)") +
  xlab("") + 
  ylab("") + 
  ggtitle("GSEA HCTbc: Resistance Genes") +
  theme_minimal() +
  theme(text = element_text(size = 7))

H_DT_fgsea_comp_plot_2 <- ggplot(data = H_DT_res_fgsea_df, 
                                 aes(x = DT3_NES, y = DT4_NES,
                                     fill = DT3_NES)) +
  geom_point(shape = 21, colour = "black", size = 7, alpha = 0.6) +
  geom_text(aes(label = gsub("_", "\n", H_DT_res_fgsea_df$pathway)),
            size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_distiller(palette = "Spectral",
                     direction = -1,
                     name = "DT3 NES") + 
  xlab("\nDT3: NES") + 
  ylab("DT4: NES\n") + 
  theme_minimal()


H_DT_fgsea_comp_plot <- cowplot::plot_grid(H_DT_fgsea_comp_plot_1,
                                           H_DT_fgsea_comp_plot_2,
                                           axis = c("lb"),
                                           align = "h",
                                           rel_widths = c(0.4, 0.5))


plot_dir <- "path/to/Barcode_Analysis/QR_10X_scRNA/DE_Analysis_Plots"

if(dir.exists(plot_dir)==F){
  dir.create(plot_dir)
  setwd(plot_dir)
}else{
  setwd(plot_dir)
}

p_height <- 8
p_width <- 15

ggsave2("QR_HCTbc_res_fgsea_cor_comp_plot.pdf",
        plot = H_DT_fgsea_comp_plot,
        height = p_height,
        width = p_width,
        device = cairo_pdf)

ggsave2("QR_HCTbc_res_fgsea_cor_comp_plot.jpg",
        plot = H_DT_fgsea_comp_plot,
        height = p_height,
        width = p_width)


# HCTbc DS FGSEA Comparison
###########################

H_DT3_per_fgsea_df <- H_DT3_per_fgsea %>% 
  rename_with(~ if_else(.x != "pathway", paste0("DT3_", .x), .x))

H_DT4_per_fgsea_df <- H_DT4_per_fgsea %>% 
  rename_with(~ if_else(.x != "pathway", paste0("DT4_", .x), .x))

H_DT_per_fgsea_df <- dplyr::full_join(H_DT3_per_fgsea_df, H_DT4_per_fgsea_df) %>%
  subset(DT3_padj <= 0.05 | DT4_padj <= 0.05) %>% 
  dplyr::mutate(pathway = gsub("HALLMARK_", "", pathway))

H_DS_fgsea_comp_plot_1 <- ggplot() + 
  geom_point(data = H_DT_per_fgsea_df, aes(x = "QR_HCTbc_DS\n(vs DT3)",
                                           y = pathway,
                                           fill = DT3_NES,
                                           size = -log(DT3_padj)),
             shape = 21, colour = "black", alpha = 0.8) + 
  geom_point(data = H_DT_per_fgsea_df, aes(x = "QR_HCTbc_DS\n(vs DT4)",
                                         y = pathway,
                                         fill = DT4_NES,
                                         size = -log(DT4_padj)),
             shape = 21, colour = "black", alpha = 0.8) + 
  scale_fill_distiller(palette = "Spectral",
                       direction = -1,
                       name = "NES") + 
  scale_size_continuous(range = c(1.0, 12.0),
                        limits = c(0.0, 70.0),
                        name = "-log(p-value, adjusted)") +
  xlab("") + 
  ylab("") + 
  ggtitle("GSEA HCTbc: Persister Genes") +
  theme_minimal()

H_DS_fgsea_comp_plot_2 <- ggplot(data = H_DT_per_fgsea_df, 
                                 aes(x = DT3_NES, y = DT4_NES,
                                     fill = DT3_NES)) +
  geom_point(shape = 21, colour = "black", size = 7, alpha = 0.6) +
  geom_text(aes(label = gsub("_", "\n", H_DT_per_fgsea_df$pathway)),
            size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_distiller(palette = "Spectral",
                     direction = -1,
                     name = "DS (vs DT3)\nNES") + 
  xlab("\nDS (vs DT3): NES") + 
  ylab("DS (vs DT4): NES\n") + 
  theme_minimal()


H_DS_fgsea_comp_plot <- cowplot::plot_grid(H_DS_fgsea_comp_plot_1,
                                           H_DS_fgsea_comp_plot_2,
                                           axis = c("lb"),
                                           align = "h",
                                           rel_widths = c(0.5, 0.5))

plot_dir <- "path/to/Barcode_Analysis/QR_10X_scRNA/DE_Analysis_Plots"

if(dir.exists(plot_dir)==F){
  dir.create(plot_dir)
  setwd(plot_dir)
}else{
  setwd(plot_dir)
}

p_height <- 8
p_width <- 15

ggsave2("QR_HCTbc_per_fgsea_cor_comp_plot.pdf",
        plot = H_DS_fgsea_comp_plot,
        height = p_height,
        width = p_width,
        device = cairo_pdf)

ggsave2("QR_HCTbc_per_fgsea_cor_comp_plot.jpg",
        plot = H_DS_fgsea_comp_plot,
        height = p_height,
        width = p_width)

#######
# SW6bc
#######

# DT1
#####

S_DT1_res_plot <- ggplot() + 
  geom_point(data = subset(S_DT1_res_df, comp_both_filter==F), 
                           aes(x = comp_1_logFC, y = comp_2_logFC,
                               colour = "F"),
             size = 3, alpha = 0.4) + 
  geom_point(data = subset(S_DT1_res_df, comp_both_filter==T), 
                           aes(x = comp_1_logFC, y = comp_2_logFC,
                               colour = "T"),
             size = 3, alpha= 0.8) + 
  xlab("logFC(DT1 vs POT)") + 
  ylab("logFC(DS vs POT)") + 
  scale_colour_manual(values = c("grey", Sdt_pal[1])) +
  theme_minimal() + 
  theme(text = element_text(size = 18)) +
  ggtitle("Resistance Treatment Response") 

S_DT1_acu_plot <- ggplot() + 
  geom_point(data = subset(S_DT1_acu_df, comp_both_filter==F), 
                           aes(x = comp_1_logFC, y = comp_2_logFC,
                               colour = "F"),
             size = 3, alpha = 0.4) + 
  geom_point(data = subset(S_DT1_acu_df, comp_both_filter==T), 
                           aes(x = comp_1_logFC, y = comp_2_logFC,
                               colour = "T"),
             size = 3, alpha= 0.8) + 
  xlab("logFC(DT1 vs POT)") + 
  ylab("logFC(DS vs POT)") + 
  scale_colour_manual(values = c("grey", "purple2")) +
  theme_minimal() + 
  theme(text = element_text(size = 18)) +
  ggtitle("Acute Treatment Response")

S_DT1_per_plot <- ggplot() + 
  geom_point(data = subset(S_DT1_per_df, comp_both_filter==F), 
                           aes(x = comp_2_logFC, y = comp_1_logFC,
                               colour = "F"),
             size = 3, alpha = 0.4) + 
  geom_point(data = subset(S_DT1_per_df, comp_both_filter==T), 
                           aes(x = comp_2_logFC, y = comp_1_logFC,
                               colour = "T"),
             size = 3, alpha = 0.8) + 
  xlab("logFC(DT1 vs POT)") + 
  ylab("logFC(DS vs POT)") + 
  scale_colour_manual(values = c("grey", Sds_pal[1])) +
  theme_minimal() + 
  theme(text = element_text(size = 18)) +
  ggtitle("Sensitive Treatment Response")

p1 <- cowplot::plot_grid(S_DT1_res_plot, 
                         S_DT1_acu_plot, 
                         S_DT1_per_plot,
                         ncol = 1)

p2 <- ggplot(data = subset(S_DT_comp_ngenes_df, samp == "S_DT1"),
       aes(x = interaction(dir, comp), 
           y = n_genes,
           fill = interaction(samp, comp))) + 
  geom_col(alpha = 0.8, colour = "black") + 
  scale_fill_manual(values = c("purple2",
                               Sds_pal[1],
                               Sdt_pal[1]),
                    name = "") +
  scale_y_continuous(limits = c(-2000, 2000)) + 
  xlab("") +
  ylab("\n\nDifferentially\nExpressed genes\n") +
  theme_minimal() +
  theme(text = element_text(size = 18))


S_DT1_comp_plot <- cowplot::plot_grid(p1, p2, nrow = 1,
                                      rel_widths = c(1.0, 0.9),
                                      axis = "tb", align = "h")


plot_dir <- "path/to/Barcode_Analysis/QR_10X_scRNA/DE_Analysis_Plots"

if(dir.exists(plot_dir)==F){
  dir.create(plot_dir)
  setwd(plot_dir)
}else{
  setwd(plot_dir)
}

p_height <- 12
p_width <- 14

ggsave2("QR_SW6bc_DT1_comp_plot.pdf",
        plot = S_DT1_comp_plot,
        height = p_height,
        width = p_width,
        device = cairo_pdf)

ggsave2("QR_SW6bc_DT1_comp_plot.jpg",
        plot = S_DT1_comp_plot,
        height = p_height,
        width = p_width)


#######
# SW6bc
#######

# DT3
#####

S_DT3_res_plot <- ggplot() + 
  geom_point(data = subset(S_DT3_res_df, comp_both_filter==F), 
                           aes(x = comp_1_logFC, y = comp_2_logFC,
                               colour = "F"),
             size = 3, alpha = 0.4) + 
  geom_point(data = subset(S_DT3_res_df, comp_both_filter==T), 
                           aes(x = comp_1_logFC, y = comp_2_logFC,
                               colour = "T"),
             size = 3, alpha= 0.8) + 
  xlab("logFC(DT3 vs POT)") + 
  ylab("logFC(DS vs POT)") + 
  scale_colour_manual(values = c("grey", Sdt_pal[4])) +
  theme_minimal() + 
  theme(text = element_text(size = 18)) +
  ggtitle("Resistance Treatment Response") 

S_DT3_acu_plot <- ggplot() + 
  geom_point(data = subset(S_DT3_acu_df, comp_both_filter==F), 
                           aes(x = comp_1_logFC, y = comp_2_logFC,
                               colour = "F"),
             size = 3, alpha = 0.4) + 
  geom_point(data = subset(S_DT3_acu_df, comp_both_filter==T), 
                           aes(x = comp_1_logFC, y = comp_2_logFC,
                               colour = "T"),
             size = 3, alpha= 0.8) + 
  xlab("logFC(DT3 vs POT)") + 
  ylab("logFC(DS vs POT)") + 
  scale_colour_manual(values = c("grey", "purple4")) +
  theme_minimal() + 
  theme(text = element_text(size = 18)) +
  ggtitle("Acute Treatment Response")

S_DT3_per_plot <- ggplot() + 
  geom_point(data = subset(S_DT3_per_df, comp_both_filter==F), 
                           aes(x = comp_2_logFC, y = comp_1_logFC,
                               colour = "F"),
             size = 3, alpha = 0.4) + 
  geom_point(data = subset(S_DT3_per_df, comp_both_filter==T), 
                           aes(x = comp_2_logFC, y = comp_1_logFC,
                               colour = "T"),
             size = 3, alpha = 0.8) + 
  xlab("logFC(DT3 vs POT)") + 
  ylab("logFC(DS vs POT)") + 
  scale_colour_manual(values = c("grey", Sds_pal[3])) +
  theme_minimal() + 
  theme(text = element_text(size = 18)) +
  ggtitle("Sensitive Treatment Response")

p1 <- cowplot::plot_grid(S_DT3_res_plot, 
                         S_DT3_acu_plot, 
                         S_DT3_per_plot,
                         ncol = 1)

p2 <- ggplot(data = subset(S_DT_comp_ngenes_df, samp == "S_DT3"),
       aes(x = interaction(dir, comp), 
           y = n_genes,
           fill = interaction(samp, comp))) + 
  geom_col(alpha = 0.8, colour = "black") + 
  scale_fill_manual(values = c("purple4",
                               Sds_pal[3],
                               Sdt_pal[4])) +
  scale_y_continuous(limits = c(-2000, 2000)) + 
  xlab("") +
  ylab("Differentially\nExpressed genes\n") +
  theme_minimal() +
  theme(text = element_text(size = 18))


S_DT3_comp_plot <- cowplot::plot_grid(p1, p2, nrow = 1,
                                      rel_widths = c(1.0, 0.9),
                                      axis = "tb", align = "h")


plot_dir <- "path/to/Barcode_Analysis/QR_10X_scRNA/DE_Analysis_Plots"

if(dir.exists(plot_dir)==F){
  dir.create(plot_dir)
  setwd(plot_dir)
}else{
  setwd(plot_dir)
}

p_height <- 12
p_width <- 14

ggsave2("QR_SW6bc_DT3_comp_plot.pdf",
        plot = S_DT3_comp_plot,
        height = p_height,
        width = p_width,
        device = cairo_pdf)

ggsave2("QR_SW6bc_DT3_comp_plot.jpg",
        plot = S_DT3_comp_plot,
        height = p_height,
        width = p_width)

#######
# HCTbc
#######

# DT3
#####

H_DT3_res_plot <- ggplot() + 
  geom_point(data = subset(H_DT3_res_df, comp_both_filter==F), 
                           aes(x = comp_1_logFC, y = comp_2_logFC,
                               colour = "F"),
             size = 3, alpha = 0.4) + 
  geom_point(data = subset(H_DT3_res_df, comp_both_filter==T), 
                           aes(x = comp_1_logFC, y = comp_2_logFC,
                               colour = "T"),
             size = 3, alpha= 0.8) + 
  xlab("logFC(DT3 vs POT)") + 
  ylab("logFC(DS vs POT)") + 
  scale_colour_manual(values = c("grey", Hdt_pal[1])) +
  theme_minimal() + 
  theme(text = element_text(size = 18)) +
  ggtitle("Resistance Treatment Response") 

H_DT3_acu_plot <- ggplot() + 
  geom_point(data = subset(H_DT3_acu_df, comp_both_filter==F), 
                           aes(x = comp_1_logFC, y = comp_2_logFC,
                               colour = "F"),
             size = 3, alpha = 0.4) + 
  geom_point(data = subset(H_DT3_acu_df, comp_both_filter==T), 
                           aes(x = comp_1_logFC, y = comp_2_logFC,
                               colour = "T"),
             size = 3, alpha= 0.8) + 
  xlab("logFC(DT3 vs POT)") + 
  ylab("logFC(DS vs POT)") + 
  scale_colour_manual(values = c("grey", "indianred1")) +
  theme_minimal() +
  theme(text = element_text(size = 18)) +
  ggtitle("Acute Treatment Response")

H_DT3_per_plot <- ggplot() + 
  geom_point(data = subset(H_DT3_per_df, comp_both_filter==F), 
                           aes(x = comp_2_logFC, y = comp_1_logFC,
                               colour = "F"),
             size = 3, alpha = 0.4) + 
  geom_point(data = subset(H_DT3_per_df, comp_both_filter==T), 
                           aes(x = comp_2_logFC, y = comp_1_logFC,
                               colour = "T"),
             size = 3, alpha = 0.8) + 
  xlab("logFC(DT3 vs POT)") + 
  ylab("logFC(DS vs POT)") + 
  scale_colour_manual(values = c("grey", Hds_pal[1])) +
  theme_minimal() + 
  theme(text = element_text(size = 18)) +
  ggtitle("Sensitive Treatment Response")

p1 <- cowplot::plot_grid(H_DT3_res_plot, 
                         H_DT3_acu_plot, 
                         H_DT3_per_plot,
                         ncol = 1)

p2 <- ggplot(data = subset(H_DT_comp_ngenes_df, samp == "H_DT3"),
       aes(x = interaction(dir, comp), 
           y = n_genes,
           fill = interaction(samp, comp))) + 
  geom_col(alpha = 0.8, colour = "black") + 
  scale_fill_manual(values = c("indianred1",
                               Hds_pal[1],
                               Hdt_pal[1])) +
  scale_y_continuous(limits = c(-1500, 1500)) + 
  xlab("") +
  ylab("Differentially\nExpressed genes\n") +
  theme_minimal() +
  theme(text = element_text(size = 18))


H_DT3_comp_plot <- cowplot::plot_grid(p1, p2, nrow = 1,
                                      rel_widths = c(1.0, 0.9),
                                      axis = "tb", align = "h")


plot_dir <- "path/to/Barcode_Analysis/QR_10X_scRNA/DE_Analysis_Plots"

if(dir.exists(plot_dir)==F){
  dir.create(plot_dir)
  setwd(plot_dir)
}else{
  setwd(plot_dir)
}

p_height <- 12
p_width <- 14

ggsave2("QR_HCTbc_DT3_comp_plot.pdf",
        plot = H_DT3_comp_plot,
        height = p_height,
        width = p_width,
        device = cairo_pdf)

ggsave2("QR_HCTbc_DT3_comp_plot.jpg",
        plot = H_DT3_comp_plot,
        height = p_height,
        width = p_width)

#######
# HCTbc
#######

# DT4
#####

H_DT4_res_plot <- ggplot() + 
  geom_point(data = subset(H_DT4_res_df, comp_both_filter==F), 
                           aes(x = comp_1_logFC, y = comp_2_logFC,
                               colour = "F"),
             size = 3, alpha = 0.4) + 
  geom_point(data = subset(H_DT4_res_df, comp_both_filter==T), 
                           aes(x = comp_1_logFC, y = comp_2_logFC,
                               colour = "T"),
             size = 3, alpha= 0.8) + 
  xlab("logFC(DT4 vs POT)") + 
  ylab("logFC(DS vs POT)") + 
  scale_colour_manual(values = c("grey", Hdt_pal[4])) +
  theme_minimal() + 
  theme(text = element_text(size = 18)) +
  ggtitle("Resistance Treatment Response") 

H_DT4_acu_plot <- ggplot() + 
  geom_point(data = subset(H_DT4_acu_df, comp_both_filter==F), 
                           aes(x = comp_1_logFC, y = comp_2_logFC,
                               colour = "F"),
             size = 3, alpha = 0.4) + 
  geom_point(data = subset(H_DT4_acu_df, comp_both_filter==T), 
                           aes(x = comp_1_logFC, y = comp_2_logFC,
                               colour = "T"),
             size = 3, alpha= 0.8) + 
  xlab("logFC(DT4 vs POT)") + 
  ylab("logFC(DS vs POT)") + 
  scale_colour_manual(values = c("grey", "indianred3")) +
  theme_minimal() +
  theme(text = element_text(size = 18)) +
  ggtitle("Acute Treatment Response")

H_DT4_per_plot <- ggplot() + 
  geom_point(data = subset(H_DT4_per_df, comp_both_filter==F), 
                           aes(x = comp_2_logFC, y = comp_1_logFC,
                               colour = "F"),
             size = 3, alpha = 0.4) + 
  geom_point(data = subset(H_DT4_per_df, comp_both_filter==T), 
                           aes(x = comp_2_logFC, y = comp_1_logFC,
                               colour = "T"),
             size = 3, alpha = 0.8) + 
  xlab("logFC(DT4 vs POT)") + 
  ylab("logFC(DS vs POT)") + 
  scale_colour_manual(values = c("grey", Hds_pal[3])) +
  theme_minimal() + 
  theme(text = element_text(size = 18)) +
  ggtitle("Sensitive Treatment Response")

p1 <- cowplot::plot_grid(H_DT4_res_plot, 
                         H_DT4_acu_plot, 
                         H_DT4_per_plot,
                         ncol = 1)

p2 <- ggplot(data = subset(H_DT_comp_ngenes_df, samp == "H_DT4"),
       aes(x = interaction(dir, comp), 
           y = n_genes,
           fill = interaction(samp, comp))) + 
  geom_col(alpha = 0.8, colour = "black") + 
  scale_fill_manual(values = c("indianred3",
                               Hds_pal[3],
                               Hdt_pal[4])) +
  scale_y_continuous(limits = c(-1500, 1500)) + 
  xlab("") +
  ylab("Differentially\nExpressed genes\n") +
  theme_minimal() + 
  theme(text = element_text(size = 18))


H_DT4_comp_plot <- cowplot::plot_grid(p1, p2, nrow = 1,
                                      rel_widths = c(1.0, 0.9),
                                      axis = "tb", align = "h")


plot_dir <- "path/to/Barcode_Analysis/QR_10X_scRNA/DE_Analysis_Plots"

if(dir.exists(plot_dir)==F){
  dir.create(plot_dir)
  setwd(plot_dir)
}else{
  setwd(plot_dir)
}

p_height <- 12
p_width <- 14

ggsave2("QR_HCTbc_DT4_comp_plot.pdf",
        plot = H_DT4_comp_plot,
        height = p_height,
        width = p_width,
        device = cairo_pdf)

ggsave2("QR_HCTbc_DT4_comp_plot.jpg",
        plot = H_DT4_comp_plot,
        height = p_height,
        width = p_width)


####################################
# DT Shared vs Unique Venn Diagrams:
####################################

plot_dir <- "path/to/Barcode_Analysis/QR_10X_scRNA/DE_Analysis_Plots"

if(dir.exists(plot_dir)==F){
  dir.create(plot_dir)
  setwd(plot_dir)
}else{
  setwd(plot_dir)
}

#######
# SW6bc
#######

S_gene_overlaps <- c(DT1 = setdiff(S_DT1_res_genes,
                                   S_DT3_res_genes) %>% length(),
                     DT3 = setdiff(S_DT3_res_genes,
                                   S_DT1_res_genes) %>% length(),
                     "DT1&DT3" = intersect(S_DT1_res_genes,
                                           S_DT3_res_genes) %>% length())

S_venn_plot <- plot(euler(S_gene_overlaps), 
                     shape = "ellipse",
                     quantities = T,
                     fills = c(Sdt_pal[c(1, 4)]),
                     alpha = 0.6,
                     lwd = 2.0)

pdf(file = "QR_SW6bc_DT_venn_plot.pdf",
    height = 6, 
    width = 10)


dev.off()

#######
# HCTbc
#######

H_gene_overlaps <- c(DT3 = setdiff(H_DT3_res_genes,
                                   H_DT4_res_genes) %>% length(),
                     DT4 = setdiff(H_DT4_res_genes,
                                   H_DT3_res_genes) %>% length(),
                     "DT3&DT4" = intersect(H_DT3_res_genes,
                                           H_DT4_res_genes) %>% length())

H_venn_plot <- plot(euler(H_gene_overlaps), 
                     shape = "ellipse",
                     quantities = T,
                     fills = c(Hdt_pal[c(1, 4)]),
                     alpha = 0.6,
                     lwd = 2.0)

pdf(file = "QR_HCTbc_DT_venn_plot.pdf",
    height = 6, 
    width = 10)


dev.off()
