# IC50 Analysis Script

# Clear the workspace
rm(list = ls())

################################################################################
# Load required libraries
################################################################################

library(drc)
library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

################################################################################
# Colour palettes
################################################################################

Hco_pal <- colorRampPalette(brewer.pal("YlGnBu", n = 9)[3:6])(4)[c(3, 3, 3, 3)]
Hdt_pal <- colorRampPalette(brewer.pal("RdPu", n = 9)[6:9])(4)[c(3, 3, 3, 3)]
Sco_pal <- colorRampPalette(brewer.pal("PuBu", n = 9)[3:6])(4)[c(3, 3, 3, 3)]
Sdt_pal <- colorRampPalette(brewer.pal("OrRd", n = 9)[6:9])(4)[c(3, 3, 3, 3)]

Hds_pal <- colorRampPalette(brewer.pal("YlOrBr", n = 9)[3:6])(4)[c(3, 3, 3, 3)]
Sds_pal <- colorRampPalette(brewer.pal("YlOrRd", n = 9)[3:6])(4)[c(3, 3, 3, 3)]

HPOT_col <- colorRampPalette(brewer.pal("Greens", n = 9)[4])(1)
SPOT_col <- colorRampPalette(brewer.pal("Greens", n = 9)[7])(1)

################################################################################
# Functions
################################################################################

#' Get IC50 DataFrame
#' 
#' Reads and processes a CSV file to extract IC50-related data.
#' 
#' @param csv_name Path to the CSV file.
#' @param dat_type Type of data to return ("csv", "newdat", "EDs", "IC50").
#' @param dat_format Format of the data (1 or 2).
#' @return DataFrame based on the requested `dat_type`.
get_IC50_df <- function(csv_name, dat_type, dat_format = 1) {
  drug <- "Fluorouracil"
  csv <- read.csv(csv_name, stringsAsFactors = FALSE)

  # Extract metadata
  assay_date <- sub("(\\d{2}_\\d{2}_\\d{2}).*", "\\1", csv_name)
  samp_name <- sub("\\d{2}_\\d{2}_\\d{2}_(.*).csv", "\\1", csv_name)

  # Process raw data
  if (dat_format == 2) {
    csv <- csv[2:7, ]
    csv <- csv[, -c(1, 2, 13)]
    csv <- data.frame(sapply(csv, as.numeric))
    csv[, 1:9] <- csv[, 1:9] - mean(csv$X11)
  } else {
    csv <- csv[13:18, ]
    csv <- csv[, -c(1:2, 13)]
    csv <- data.frame(sapply(csv, as.numeric))
    csv[, 1:9] <- csv[, 1:9] - mean(csv$X.10)
  }

  csv[, 1:9][csv[, 1:9] < 0] <- 0.0
  csv[, 1:9] <- csv[, 1:9] / mean(csv[[9]]) * 100
  csv <- csv[, 1:9]
  names(csv) <- concs
  csv <- gather(csv, "Concn", "Response", 1:9)
  csv$Concn <- as.numeric(csv$Concn)

  # Fit the dose-response curve
  drug_curve <- drm(Response ~ Concn, data = csv, fct = LL.3())
  EDs <- data.frame(ED(drug_curve, c(50, 60, 70, 80), interval = "delta"))
  EDs["samp_name"] <- samp_name
  EDs["ED"] <- as.numeric(sub(".*:(\\d{2})", rownames(EDs), replacement = "\\1"))

  # Generate prediction data
  newdat <- expand.grid(Concn = exp(seq(log(0.01), log(2000), length = 1000)))
  p_mox <- predict(drug_curve, newdata = newdat, interval = "confidence")
  newdat$p <- p_mox[, 1]
  newdat$pmin <- p_mox[, 2]
  newdat$pmax <- p_mox[, 3]
  newdat["samp_name"] <- samp_name

  IC50 <- newdat[which.min(abs(50 - newdat$p)), "Concn"]

  if (dat_type == "csv") return(csv)
  if (dat_type == "newdat") return(newdat)
  if (dat_type == "EDs") return(EDs)
  if (dat_type == "IC50") return(data.frame(IC50 = IC50, samp_name = samp_name))
}

#' Generate Summary DataFrame
#' 
#' Aggregates mean and standard deviation of cell viability by concentration.
#' 
#' @param samp_csv DataFrame containing raw sample data.
#' @return DataFrame with summary statistics.
get_summ_df <- function(samp_csv) {
  plyr::ddply(samp_csv, .(Concn, samp_name), summarise, mean = mean(Cell_Viability), sd = sd(Cell_Viability))
}

#' Plot Combined IC50 Data
#' 
#' Creates a ggplot visualisation of the IC50 data, including confidence intervals.
#' 
#' @param samp_csv Sample raw data.
#' @param samp_newdat Sample prediction data.
#' @param samp_IC50 Sample IC50 values.
#' @param samp_pal Colour palette for the sample.
#' @param POT_csv POT sample raw data.
#' @param POT_newdat POT sample prediction data.
#' @param POT_IC50 POT sample IC50 values.
#' @param POT_pal Colour palette for the POT sample.
#' @return ggplot object.
plot_comb_IC50s <- function(samp_csv, samp_newdat, samp_IC50, samp_pal, 
                            POT_csv, POT_newdat, POT_IC50, POT_pal) {
  samp_summ <- get_summ_df(samp_csv)
  POT_summ <- get_summ_df(POT_csv)

  # Add the POT summary and newdat to the sample dataframe
  POT_summ <- POT_summ[, c("Concn", "mean", "sd")]
  colnames(POT_summ) <- c("Concn", "POT_mean", "POT_sd")
  POT_newdat <- POT_newdat[, c("Concn", "p")]
  colnames(POT_newdat) <- c("Concn", "POT_p")

  samp_summ <- left_join(samp_summ, POT_summ, by = "Concn")
  samp_newdat <- left_join(samp_newdat, POT_newdat, by = "Concn")

  ggplot(data = samp_summ, aes(x = Concn + 1e-03, y = mean, colour = samp_name, fill = samp_name)) +
    geom_line(data = samp_newdat, aes(x = Concn + 1e-03, y = p), size = 1) +
    geom_errorbar(aes(x = Concn + 1e-03, ymin = mean - sd, ymax = mean + sd), width = 0.05) +
    geom_point(shape = 21, colour = "black", size = 3) +
    geom_line(data = samp_newdat, aes(x = Concn + 1e-03, y = POT_p), size = 1, colour = POT_pal, alpha = 0.7) +
    geom_errorbar(aes(x = Concn + 1e-03, ymin = POT_mean - sd, ymax = POT_mean + sd), width = 0.05, colour = POT_pal, alpha = 0.7) +
    geom_point(aes(x = Concn + 1e-03, y = POT_mean), shape = 21, colour = "black", size = 3, fill = POT_pal, alpha = 0.7) +
    geom_vline(aes(xintercept = POT_IC50$IC50 + 1e-03), linetype = "dashed", size = 1, colour = "black", alpha = 0.6) +
    geom_vline(data = samp_IC50, aes(xintercept = IC50 + 1e-03), linetype = "dashed", size = 1, colour = "red") +
    scale_x_log10(limits = c(1, 300),
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    annotation_logticks(base = 10, sides = "b") +
    scale_y_continuous(limits = c(0, 120), breaks = seq(0, 100, 25), labels = seq(0, 100, 25)) +
    facet_wrap(~samp_name, ncol = 1) +
    scale_colour_manual(values = samp_pal, name = "Sample") +
    scale_fill_manual(values = samp_pal, name = "Sample") +
    theme_minimal() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none") +
    xlab("\nlog10[5-Fu]") +
    ylab("Viability (%)\n")
}

################################################################################
# Main Script
################################################################################

# Define input files and concentrations
input_files <- list.files("/path/to/csv/files", pattern = "*.csv", full.names = TRUE)
concs <- c(200, 100, 50, 20, 10, 5, 2, 1, 0)

# Example usage
example_file <- input_files[1]
example_csv <- get_IC50_df(example_file, dat_type = "csv", dat_format = 1)
example_newdat <- get_IC50_df(example_file, dat_type = "newdat", dat_format = 1)
example_IC50 <- get_IC50_df(example_file, dat_type = "IC50", dat_format = 1)
example_plot <- plot_comb_IC50s(example_csv, example_newdat, example_IC50, Hco_pal, 
                                example_csv, example_newdat, example_IC50, HPOT_col)

# Save the plot
ggsave("example_plot.pdf", example_plot, height = 6, width = 4)
