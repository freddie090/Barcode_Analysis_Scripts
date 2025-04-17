# QR Experiment - Single Colony Expansion Analysis
# Deriving growth rates via the time-series confluence data from the Incucyte analysis.

rm(list = ls())

# Load necessary libraries
library(tidyr)
library(plyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(RColorBrewer)

# Define a generic working directory
set_working_directory <- function(path = "./data") {
  if (!dir.exists(path)) {
    stop("Directory does not exist: ", path)
  }
  setwd(path)
}

# Load data files matching a pattern
load_data_files <- function(pattern) {
  file_names <- grep(pattern, list.files(), value = TRUE)
  if (length(file_names) == 0) {
    stop("No files matching pattern found: ", pattern)
  }
  lapply(file_names, read.table, skip = 1, header = TRUE)
}

# Extract sample names from file names
extract_sample_names <- function(file_names) {
  sub("(QR_.*_P\\d{1})_sc.*", "\\1", file_names)
}

# Define selected wells for each sample
create_well_data <- function() {
  list(
    data.frame(samp = "QR_HCTbc_CO3_P4", well = c("B2", "B4", "B11", "C6", "C9", "D2", "D10", "E8", "F5", "F10", "G6", "G9")),
    data.frame(samp = "QR_HCTbc_CO4_P4", well = c("B2", "B3", "B7", "B11", "C2", "C5", "D4", "E7", "E8", "F3", "F7", "G2")),
    data.frame(samp = "QR_HCTbc_DT3_P4", well = c("B1", "B3", "B4", "B8", "B10", "C3", "C4", "C5", "C7", "C9", "D4", "D6")),
    data.frame(samp = "QR_HCTbc_DT4_P4", well = c("B3", "B5", "C2", "C9", "D5", "D7", "D11", "E2", "E6", "E7", "F7", "G7"))
  ) %>% bind_rows()
}

# Process loaded data
process_data <- function(data_list, sample_names, well_data) {
  for (i in seq_along(sample_names)) {
    data_list[[i]] <- data_list[[i]][, !grepl("Date|Time", colnames(data_list[[i]]))]
    data_list[[i]] <- subset(data_list[[i]], Elapsed <= 320)
    data_list[[i]] <- gather(data_list[[i]], key = "well", value = "confluence", -Elapsed)
    data_list[[i]]["samp"] <- sample_names[[i]]
    data_list[[i]] <- data_list[[i]][data_list[[i]]$well %in% subset(well_data, samp == sample_names[[i]])$well, ]
  }
  bind_rows(data_list)
}

# Adjust elapsed time and save data
adjust_and_save_data <- function(data, output_path = "processed_confluence_data.csv") {
  data$Elapsed <- data$Elapsed + 144
  write.csv(data, output_path, row.names = FALSE)
  data
}

# Fit linear models for growth rate estimation
fit_growth_models <- function(data, separate_intercepts = TRUE) {
  samp <- unique(data$samp)
  if (length(samp) > 1) stop("There should only be one unique sample in this dataframe.")
  
  if (separate_intercepts) {
    mod_fit_df <- plyr::ddply(data, .(well), function(x) {
      t <- x$Elapsed
      y <- x$confluence + 1e-02
      log_y <- log(y)
      mod_fit <- lm(log_y ~ t)
      data.frame(intercept = coef(mod_fit)[[1]], 
                 gradient = coef(mod_fit)[[2]])
    })
  } else {
    dat_df <- data.frame(t = data$Elapsed, y = data$confluence, well = data$well)
    dat_df$y <- dat_df$y + 1e-02
    dat_df$log_y <- log(dat_df$y)
    mod_fit <- lm(log_y ~ t:well, data = dat_df)
    fit_inter <- coef(mod_fit)[[1]]
    mod_df <- data.frame(gradient = coef(mod_fit)[-1])
    mod_df$well <- sub(".*well(\\D{1}\\d{1,2})", "\\1", rownames(mod_df))
    mod_df$intercept <- fit_inter
    rownames(mod_df) <- NULL
    mod_fit_df <- mod_df
  }
  return(mod_fit_df)
}

# Generate simulated data based on model fits
simulate_growth_data <- function(models, time_seq) {
  plyr::ddply(models, .(samp, well), function(row) {
    intercept <- row$intercept
    gradient <- row$gradient
    log_y <- intercept + gradient * time_seq
    data.frame(t = time_seq, log_y = log_y, y = exp(log_y))
  })
}

# Save growth model fits to CSV
save_growth_models <- function(models, output_path) {
  write.csv(models, output_path, row.names = FALSE)
}


# Plot and save confluence vs time for a specific sample
well_confluence_plot <- function(comb_out_df, comb_sim_df, samp_name, log_y, output_dir = "./plots") {
  out_samp_df <- subset(comb_out_df, samp == samp_name)
  sim_samp_df <- subset(comb_sim_df, samp == samp_name)
  pal12 <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928')
  
  p <- ggplot(data = out_samp_df, aes(x = Elapsed, y = confluence, fill = well)) + 
    geom_point(shape = 21, colour = "black", alpha = 0.8, size = 2) + 
    geom_line(data = sim_samp_df, aes(x = t, y = y, colour = well), size = 1) + 
    facet_grid(~well) + 
    theme_minimal() + 
    theme(legend.position = "none") +
    ggtitle(samp_name) + 
    scale_y_continuous(limits = c(1e-02, 20)) +
    scale_color_manual(values = pal12) + 
    scale_fill_manual(values = pal12) + 
    theme(axis.text.x = element_blank(),
          text = element_text(size = 22)) + 
    xlab("Elapsed (hrs)") + 
    ylab("Confluence (%)")
  
  if (log_y) {
    p <- p + scale_y_log10(limits = c(1e-03, 20))
  }
  
  ggsave(file.path(output_dir, paste0("fit_comp_", samp_name, ".pdf")), p, height = 8, width = 10)
  print(p)
}


# Plot gradient comparison
plot_gradient_comparison <- function(models, output_path = "./plots/gradient_comparison.pdf") {
  Hco_pal <- colorRampPalette(brewer.pal("YlGnBu", n = 9)[3:6])(4)
  Hdt_pal <- colorRampPalette(brewer.pal("RdPu", n = 9)[6:9])(4)

  gradient_plot <- ggplot(models, aes(x = samp, y = gradient, fill = samp)) +
    geom_boxplot(alpha = 0.5) +
    geom_jitter(shape = 21, colour = "black", size = 2, alpha = 0.7, width = 0.1) +
    theme_minimal() +
    scale_y_continuous(limits = c(0.0000, 0.035)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = c(Hco_pal[c(1,4)], Hdt_pal[c(1,4)]), name = "Sample") +
    xlab("\nSample") +
    ylab("Gradient Fit\n") +
    stat_compare_means(
      comparisons = list(
        c("HCTbc CO3 P4", "HCTbc CO4 P4"),
        c("HCTbc CO3 P4", "HCTbc DT3 P4"),
        c("HCTbc CO3 P4", "HCTbc DT4 P4"),
        c("HCTbc CO4 P4", "HCTbc DT3 P4"),
        c("HCTbc CO4 P4", "HCTbc DT4 P4"),
        c("HCTbc DT3 P4", "HCTbc DT4 P4")
      ),
      method = "wilcox.test",
      correction = "BH",
      label = "p.signif"
    )
  
  ggsave(output_path, gradient_plot, height = 3, width = 8)
  print(gradient_plot)
}

# Main workflow
main <- function() {
  #set_working_directory("./data")
  set_working_directory("C:/Users/fwhiting/My Drive/Lab/Experiments/QR_Single_Cell_Sorts/Incucyte_Analysis/")
  file_pattern <- "QR_.*_Analysis1.txt"
  data_list <- load_data_files(file_pattern)
  sample_names <- extract_sample_names(grep("QR_.*_Analysis1.txt", 
                                            list.files(), value = T))
  well_data <- create_well_data()
  processed_data <- process_data(data_list, sample_names, well_data)
  adjusted_data <- adjust_and_save_data(processed_data)

  growth_models <- plyr::ddply(adjusted_data, .(samp), fit_growth_models,
                               separate_intercepts = F)
  save_growth_models(growth_models, "growth_models.csv")

  
  time_seq <- seq(0, max(adjusted_data$Elapsed), by = 2)
  simulated_data <- simulate_growth_data(growth_models, time_seq)

  # Generate and save individual sample plots
  unique_samples <- unique(adjusted_data$samp)
  for (samp in unique_samples) {
    well_confluence_plot(adjusted_data, simulated_data, samp, log_y = T)
  }
  
  plot_gradient_comparison(growth_models, output_path = "./plots/gradient_comparison.pdf")
  
}

# Run the script
main()
