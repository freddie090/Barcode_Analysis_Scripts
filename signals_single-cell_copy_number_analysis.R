# Clear environment
rm(list = ls())

# Load required libraries
library(tidyr)
library(plyr)
library(dplyr)
library(ggplot2)
library(signals)
library(RColorBrewer)
library(cowplot)

# Set working directory
setwd("/path/to/Barcode_Analysis/QR_DLPs/")

################################################################################

#' Create a Minimal Custom ggplot2 Theme
#'
#' @return A ggplot2 theme object.
theme_BARCODE <- function() {
  theme_minimal() +
    theme(
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      text = element_text(size = 18)
    )
}

################################################################################

# Define colour palettes for plotting
Hco_pal <- colorRampPalette(brewer.pal("YlGnBu", n = 9)[3:6])(4)
Hdt_pal <- colorRampPalette(brewer.pal("RdPu", n = 9)[6:9])(4)
Sco_pal <- colorRampPalette(brewer.pal("PuBu", n = 9)[3:6])(4)
Sdt_pal <- colorRampPalette(brewer.pal("OrRd", n = 9)[6:9])(4)

Hds_pal <- colorRampPalette(brewer.pal("YlOrBr", n = 9)[3:5])(3)
Sds_pal <- colorRampPalette(brewer.pal("YlOrRd", n = 9)[3:5])(3)

HPOT_col <- colorRampPalette(brewer.pal("Greens", n = 9)[4])(1)
SPOT_col <- colorRampPalette(brewer.pal("Greens", n = 9)[7])(1)

HCT_samp_colvals <- c(
  "POTs_EP1" = HPOT_col,
  "CO3_P4_EP1" = Hco_pal[1],
  "CO4_P4_EP1" = Hco_pal[4],
  "DSI_P1_EP1" = Hds_pal[1],
  "DSII_P1_EP1" = Hds_pal[3],
  "DT3_P4_EP1" = Hdt_pal[1],
  "DT4_P4_EP1" = Hdt_pal[4]
)

SW6_samp_colvals <- c(
  "POTs_EP1" = SPOT_col,
  "CO1_P4_EP1" = Sco_pal[1],
  "CO3_P4_EP1" = Sco_pal[4],
  "DSI_P4_EP1" = Sds_pal[1],
  "DSII_P1_EP1" = Sds_pal[3],
  "DT1_P4_EP1" = Sdt_pal[1],
  "DT3_P4_EP1" = Sdt_pal[4]
)

################################################################################

# Load chromosome data
chrom_df <- read.csv("/path/to/Barcode_Analysis/Human_Chromosome_Arm_Lengths.csv")

# Check if already filtered:
if (file.exists("QR_SW6bc_DLP_filtered_reads.csv") == FALSE) {

  DLP_csvs <- grep("_reads.csv", list.files(), value = TRUE)
  DLP_dfs <- lapply(DLP_csvs, read.csv)

  # Merge for filtering:
  DLP_df <- bind_rows(DLP_dfs)
  rm(DLP_dfs)

  # Get read counts for each cell:
  tot_read_df <- plyr::ddply(DLP_df, .(cell_id), summarise, 
                             tot_reads = sum(reads))

  # Pick those with >1e+05 reads:
  sub_cells <- subset(tot_read_df, tot_reads >= 1e+05)$cell_id
  DLP_df <- subset(DLP_df, cell_id %in% sub_cells)
  DLP_df <- subset(DLP_df, !is.na(copy))

  # Compute sub-cell groups for chromosomes:
  avg_ppc_df <- plyr::ddply(DLP_df, .(cell_id, chr), summarise, 
                            m_copy = mean(copy),
                            m_state = mean(state))

  sub_cells_1 <- unique(subset(avg_ppc_df, chr == 3 & m_copy <= 2.6)$cell_id)
  sub_cells_2 <- unique(subset(avg_ppc_df, chr == 1 & m_copy <= 2.6)$cell_id)

  sub_cells <- sub_cells_1[(sub_cells_1 %in% sub_cells_2)]

  DLP_df <- subset(DLP_df, cell_id %in% sub_cells)

  # Change the cell_id names so they work with the signals heatmap function:
  DLP_df$cell_id <- gsub("_(?=DLP)", "-", DLP_df$cell_id, perl = TRUE)

  # Save this filtered dataframe:
  write.csv(DLP_df, "QR_SW6bc_DLP_filtered_reads.csv",
            row.names = FALSE)

} else {
  DLP_df <- read.csv("QR_SW6bc_DLP_filtered_reads.csv")
}

# Add a sample column to the dataframe:
DLP_df$sample <- sub("(QR_SW6bc_.*)-DLP_.*", replacement = "\\1", DLP_df$cell_id)

# Split out the cell_id and sample columns:
DLP_samp_split_df <- DLP_df[, c("cell_id", "sample")] %>% 
  unique() %>% 
  group_split(sample)

# Sub-sample 50 cells per sample:
DLP_samp_cells <- lapply(DLP_samp_split_df, function(x) {
  x$cell_id[sample(1:nrow(x), min(50, nrow(x)), replace = FALSE)]
}) %>% unlist()

# Now only keep the cells in the read dataframe that match these cells:
DLP_df_sub <- subset(DLP_df, cell_id %in% DLP_samp_cells)

# Split sub-sampled data by sample:
DLP_df_subs <- DLP_df_sub %>% group_by(sample) %>% group_split()

# Cluster via UMAP:
DLP_clustering <- umap_clustering(DLP_df_sub, field = "state")

# Save heatmap
pdf("QR_SW6bc_DLP_heatmap.pdf", width = 12.0, height = 15.0)
plotHeatmap(DLP_df_sub, 
            tree = DLP_clustering$tree, 
            clusters = DLP_clustering$clustering)
dev.off()

# Generate plots for consensus profile per sample:
DLP_df_subs <- DLP_df_subs[match(c("QR_SW6bc_POTs_EP1", "QR_SW6bc_CO1_P4_EP1",
                                   "QR_SW6bc_DSI_P4_EP1", "QR_SW6bc_DSII_P1_EP1",
                                   "QR_SW6bc_DT1_P4_EP1", "QR_SW6bc_DT3_P4_EP1"),
                                 lapply(DLP_df_subs, function(x) { unique(x$sample) }) %>% unlist())]

dlp_CNbins_cons <- lapply(DLP_df_subs, consensuscopynumber, cl = DLP_clustering$clustering)

dlp_CNbins_cons_plots <- lapply(dlp_CNbins_cons, plotCNprofile)

for (i in seq_along(dlp_CNbins_cons_plots)) {
  dlp_CNbins_cons_plots[[i]] <- dlp_CNbins_cons_plots[[i]] + 
    ggtitle(unique(DLP_df_subs[[i]]$sample)) + 
    theme(text = element_text(size = 18))
}

DLP_cons_plots <- cowplot::plot_grid(plotlist = dlp_CNbins_cons_plots, ncol = 2)
ggsave2("QR_SW6bc_DLP_consensus.pdf", DLP_cons_plots, width = 13.0, height = 9.0, units = "in")

# Generate UMAP by sample:
umap_results <- DLP_clustering$umapresults$embedding %>% data.frame()
umap_results$cell_id <- rownames(umap_results)
rownames(umap_results) <- NULL

umap_results$sample <- sub("QR_SW6bc_(.*)-DLP_.*", replacement = "\\1", umap_results$cell_id)

umap_plot <- ggplot(data = umap_results, aes(x = X1, y = X2, colour = sample)) + 
  geom_jitter(size = 4.0, alpha = 0.6,# colour = "black",
              height = 0.1, width = 0.1) + 
  theme_BARCODE() + 
  scale_colour_manual(values = SW6_samp_colvals,
                      name = "Sample")

ggsave2("QR_SW6bc_umap_by_samp.pdf",
       width = 10.0, height = 8.0, units = "in")

################################################################################
