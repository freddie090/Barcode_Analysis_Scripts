
################################################################################
# QR 10X scRNA Experiment: SComatic Analysis
# Freddie Whiting
# Jan 2025
################################################################################

# This script processes scRNA data from a QR Barcode experiment to identify and analyze mutations.
# It includes functions to process mutation data, generate plots of mutation counts and variant allele frequency (VAF) distributions,
# and optionally subsets the data for specific samples.
# Results are saved as visualizations to PDF files.

################################################################################

# Load required libraries
library(ggplot2)
library(dplyr)
library(cowplot)
library(RColorBrewer)
library(scales)
library(reshape2)
library(stringr)
library(tidyr)
library(plyr)
library(purrr)

# Custom Themes for Plots
#' @title Custom Barcode Plot Themes
#' @description Provides minimalistic themes for plotting.
theme_BARCODE <- function() {
  theme_minimal() +
    theme(
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

theme_BARCODE_2 <- function() {
  theme_BARCODE() +
    theme(
      panel.grid.major = element_line(size = 0.5, colour = "grey90"),
      panel.border = element_rect(colour = "black", fill = NA, size = 1.0)
    )
}

# Colour Palettes for Plotting
# Define sample colour palettes
HCT_samp_colvals <- c(
  "POTs" = brewer.pal("Greens", n = 9)[4],
  "CO3 P4" = brewer.pal("YlGnBu", n = 9)[3],
  "CO4 P4" = brewer.pal("YlGnBu", n = 9)[6],
  "DSI P1" = brewer.pal("YlOrBr", n = 9)[3],
  "DSII P1" = brewer.pal("YlOrBr", n = 9)[5],
  "DT3 P4" = brewer.pal("RdPu", n = 9)[6],
  "DT4 P4" = brewer.pal("RdPu", n = 9)[9]
)

#' @title Split by Sample
#' @description Splits a row of a data frame by sample columns and converts to long format.
#' @param DF_ROW A single row of a data frame.
#' @return A data frame with expanded rows for each sample.
split_by_sample <- function(DF_ROW) {
  SUB_DF <- apply(DF_ROW, 2, function(x) {
    strsplit(x, ",") %>% unlist()
  })
  
  if ((SUB_DF %>% lapply(length) %>% unlist() %>% max()) > 1) {
    OUT_DF <- data.frame(SUB_DF)
  } else {
    OUT_DF <- data.frame(t(SUB_DF))
  }
  
  return(OUT_DF)
}

#' @title Process Mutation Data
#' @description Processes mutation data to separate comma-separated fields and add sample count columns.
#' @param file_path Path to the mutation data file.
#' @return A list with long and wide formatted mutation data frames.
process_mutation_data <- function(file_path) {
  mut_df <- read.table(file_path, sep = "\t", header = TRUE)
  mut_df$mut_id <- 1:nrow(mut_df)
  mut_df$nmut_samps <- lapply(strsplit(mut_df$Cell_types, ","), length) %>% unlist()

  mut_dfl <- mut_df %>%
    separate_rows(
      ALT, Cell_types, Dp, Nc, Bc, Cc, VAF, CCF, BCp, CCp, Cell_type_Filter,
      sep = ",",
      convert = TRUE
    )

  mut_dfl$nmut_samps <- paste0("Found in \n", mut_dfl$nmut_samps)
  mut_dfl$Cell_types <- gsub("_EP1", "", mut_dfl$Cell_types)
  mut_dfl$Cell_types <- gsub("_", " ", mut_dfl$Cell_types)

  mut_dfw <- mut_dfl %>%
    pivot_wider(names_from = Cell_types, values_from = c("ALT", "Dp", "Nc", "Bc", "Cc", "VAF", "CCF", "BCp", "CCp"))

  list(long = mut_dfl, wide = mut_dfw)
}

#' @title Plot Mutation VAF Distribution
#' @description Plots the distribution of mutation VAFs by sample and occurrence frequency.
#' @param mut_dfl Long-format mutation data frame.
#' @param output_path Path to save the output plot.
#' @param subset_samples Optional vector of samples to subset by. Defaults to NULL.
plot_vaf_distribution <- function(mut_dfl, output_path, subset_samples = NULL) {
  if (!is.null(subset_samples)) {
    data_to_plot <- subset(mut_dfl, Cell_types %in% subset_samples & nmut_samps == "Found in \n1")
  } else {
    data_to_plot <- mut_dfl
  }

  QR_HCTbc_DT_fi1_muts_plot <- ggplot(data_to_plot, aes(x = VAF, fill = Cell_types)) +
    geom_histogram(alpha = 0.6, position = "identity", bins = 50, colour = "black") +
    facet_grid(Cell_types ~ nmut_samps) +
    theme_BARCODE_2() +
    scale_fill_manual(values = HCT_samp_colvals) +
    ylab("Mutation\nCount\n")

  ggsave(output_path, QR_HCTbc_DT_fi1_muts_plot, height = 5, width = 5)

  QR_HCTbc_DT_fi1_muts_plot
}

# Main Script Execution

# Set file path
file_path <- "QR_HCTbc.calling.step2.pass.tsv"

# Process mutation data
mutation_data <- process_mutation_data(file_path)

# Generate and save plots
plot_total_mutations(mutation_data$long, "QR_HCTbc_tot_muts_plot.pdf")

# Plot VAF Distribution without subsetting
plot_vaf_distribution(mutation_data$long, "QR_HCTbc_vaf_full_plot.pdf")

# Plot VAF Distribution with subsetting
plot_vaf_distribution(mutation_data$long, "QR_HCTbc_DT_subset_vaf_plot.pdf", subset_samples = c("DT3 P4", "DT4 P4"))
