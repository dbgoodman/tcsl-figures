#' Main script for generating figures for arrayed flow panel round 2 CAR experiments
#'
#' This script processes flow cytometry data from TCSL248 and TCSL250 experiments
#' and generates standardized plots for analysis.

# Load required libraries
library(here)
library(data.table)
library(ggplot2)
library(patchwork)

# Source the processing and plotting functions
source(here("figures", "arrayed_flowpanel_round2_cars", "data_setup", "process_flowpanel_data.R"))
source(here("figures", "arrayed_flowpanel_round2_cars", "data_setup", "plot_flowpanel.R"))
source(here("figures", "arrayed_flowpanel_round2_cars", "process_experiment.R"))
source(here("figures", "shared", "car_colors.R"))
source(here("figures", "shared", "car_sets.R"))

# Create output directories
output_dir <- here("data", "processed", "flowpanel_round2")
plots_dir_base <- here("figures", "arrayed_flowpanel_round2_cars", "plots")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir_base, recursive = TRUE, showWarnings = FALSE)

# Function to create replicate metadata for TCSL248
create_tcsl248_rep_metadata <- function() {
  rep_metadata <- data.table(
    rep = c('X', 'A', 'B', 'C'),
    rep_type = c('Unstim', '1:8', '1:8', '1:4'),
    stimmed = c(FALSE, TRUE, TRUE, TRUE)
  )
  rep_metadata[, rep_lbl := factor(rep, levels=c('X','A','B','C'),
                                  labels = c('Unstim', '1:8a', '1:8b', '1:4'))]
  return(rep_metadata)
}

# Function to create replicate metadata for TCSL250
create_tcsl250_rep_metadata <- function() {
  rep_metadata <- data.table(
    rep = c('A', 'B', 'C'),
    rep_type = c('1:4', '1:4', '1:4'),
    stimmed = c(TRUE, TRUE, TRUE)
  )
  rep_metadata[, rep_lbl := factor(rep, levels=c('A','B','C'),
                                  labels = c('1:4a', '1:4b', '1:4c'))]
  return(rep_metadata)
}

# Define bead count files
bead_count_file_tcsl248 <- here("data", "raw", "flowpanel_round2", "tcsl248", "../20241115.tcsl248_wk1_count/out_counts_clean.csv")
if (!file.exists(bead_count_file_tcsl248)) {
  message("Bead count file not found for TCSL248")
  bead_count_file_tcsl248 <- NULL
}

bead_count_file_tcsl250 <- here("data", "raw", "flowpanel_round2", "tcsl250", "bead_count.csv")
if (!file.exists(bead_count_file_tcsl250)) {
  message("Bead count file not found for TCSL250")
  bead_count_file_tcsl250 <- NULL
}

# Process TCSL248 data
tcsl248_rep_metadata <- create_tcsl248_rep_metadata()
plots_dir_tcsl248 <- file.path(plots_dir_base, "tcsl248")
df_tcsl248 <- process_experiment("tcsl248", tcsl248_rep_metadata,
                                output_dir, plots_dir_tcsl248,
                                bead_count_file_tcsl248)

# Process TCSL250 data
tcsl250_rep_metadata <- create_tcsl250_rep_metadata()
plots_dir_tcsl250 <- file.path(plots_dir_base, "tcsl250")
df_tcsl250 <- process_experiment("tcsl250", tcsl250_rep_metadata,
                                output_dir, plots_dir_tcsl250,
                                bead_count_file_tcsl250)

# Generate combined plots
message("Generating combined plots...")
plots_dir_combined <- file.path(plots_dir_base, "combined")
dir.create(plots_dir_combined, recursive = TRUE, showWarnings = FALSE)

# Add experiment identifier to each dataset
df_tcsl248[, experiment := "tcsl248"]
df_tcsl250[, experiment := "tcsl250"]

# Combine datasets
df_combined <- rbind(df_tcsl248, df_tcsl250, fill = TRUE)

# Filter out unstimulated samples for combined plots
df_combined_stim <- df_combined[stimmed == TRUE]

# # Create combined tile plots
# message("Creating combined tile plots...")
# combined_tile_plots <- create_tile_plots(df_combined_stim)
# for (name in names(combined_tile_plots)) {
#   filepath <- file.path(plots_dir_combined, paste0("combined_tile_plot_", name, ".pdf"))
#   print(filepath)
#   ggsave(filepath,
#          combined_tile_plots[[name]],
#          width = 20, height = 10, create.dir = TRUE)
# }

# # Create combined CD3 tile plots
# combined_cd3_tile_plots <- create_tile_plots(df_combined_stim, subset_filter = "CD3")
# for (name in names(combined_cd3_tile_plots)) {
#   ggsave(file.path(plots_dir_combined, paste0("combined_cd3_tile_plot_", name, ".pdf")),
#          combined_cd3_tile_plots[[name]],
#          width = 20, height = 8, create.dir = TRUE)
# }

# Function to normalize values relative to zeta CAR
normalize_to_zeta <- function(data, value_col = "freq") {
  # Create a copy to avoid modifying the original data
  norm_data <- copy(data)
  
  # Calculate zeta values for each experiment and subset
  zeta_values <- norm_data[car_name == "Zeta",
                          .(zeta_value = mean(get(value_col), na.rm = TRUE)),
                          by = .(experiment, subset)]
  
  # Merge zeta values back to the data
  norm_data <- merge(norm_data, zeta_values, by = c("experiment", "subset"))
  
  # Calculate normalized values
  norm_data[, paste0(value_col, "_norm") := get(value_col) / zeta_value]
  
  return(norm_data)
}

# Get CAR sets
car_sets <- get_car_sets()

# Define cell subsets and marker categories
cell_subsets <- c("CD3", "CD4", "CD8")
marker_categories <- list(
  exhaustion = c("PD1", "CD39", "TIM3", "LAG3"),
  differentiation = c("Naive", "Mem", "Eff", "Emra"),
  activation = c("CD25", "CD27", "CD127")
)

# Loop through each CAR set
for (car_set_name in names(car_sets)) {
  car_set <- car_sets[[car_set_name]]
  message(sprintf("Processing CAR set: %s", car_set$title_prefix))

  # Create directory for this CAR set
  car_set_dir <- file.path(plots_dir_combined, car_set_name)
  dir.create(car_set_dir, recursive = TRUE, showWarnings = FALSE)

  # Filter data for this CAR set and ensure only stimulated data is used
  df_car_set <- filter_car_set(df_combined_stim, car_set)

  # Skip if no data for this CAR set
  if (nrow(df_car_set) == 0) {
    message(sprintf("No data for CAR set: %s, skipping", car_set$title_prefix))
    next
  }

  # Define a list to store all plots for this CAR set
  plots <- list()

  # Ensure gate categories are assigned for category-based plots
  df_car_set[gate_clean %in% c('CD25','CD27','CD127'), gate_category := 'Activation']
  df_car_set[gate_clean %in% c('Eff','Emra','Mem','Naive'), gate_category := 'Differentiation']
  df_car_set[gate_clean %in% c('TIM3','LAG3','CD39','PD1'), gate_category := 'Exhaustion']

  # Create a combined CAR expression plot with normalization
  car_data <- df_car_set[subset %in% cell_subsets & gate_clean == "live_car" &
                       xor(CAR == TRUE, car_name == "U")]

  # Double-check that we're only using stimulated data
  car_data <- car_data[stimmed == TRUE]

  # Normalize to zeta
  car_data_norm <- normalize_to_zeta(car_data)

  # Create plots for each subset
  for (subset in cell_subsets) {
    subset_car_data <- car_data_norm[subset == subset]

    # Skip if no data for this subset
    if (nrow(subset_car_data) == 0) {
      message(sprintf("No CAR expression data for %s in CAR set %s, skipping",
                     subset, car_set$title_prefix))
      next
    }

    # Create the plot with normalized values and reordered x-axis
    plot_key <- paste0("car_expression_", tolower(subset), "_normalized")
    plots[[plot_key]] <- ggplot(subset_car_data) +
      geom_point(aes(x = reorder(car_name, freq_norm, mean), y = freq_norm, color = experiment, shape = rep_lbl), size = 3) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
      scale_y_continuous(name = "CAR Expression (relative to Zeta)") +
      scale_color_manual(values = c("tcsl248" = "#E69F00", "tcsl250" = "#56B4E9"),
                         labels = c("tcsl248" = "TCSL248", "tcsl250" = "TCSL250")) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            legend.position = "right") +
      labs(title = paste0(car_set$title_prefix, " - ", subset, " - CAR Expression (Normalized to Zeta)"),
           x = "CAR Construct")

    # Also create a version with mean values per CAR
    mean_car_data <- subset_car_data[, .(mean_norm = mean(freq_norm, na.rm = TRUE),
                                       se_norm = sd(freq_norm, na.rm = TRUE)/sqrt(.N)),
                                   by = .(car_name, experiment)]

    plot_key_mean <- paste0("car_expression_", tolower(subset), "_normalized_mean")
    plots[[plot_key_mean]] <- ggplot(mean_car_data) +
      geom_point(aes(x = reorder(car_name, mean_norm, mean), y = mean_norm, color = experiment), size = 3) +
      geom_errorbar(aes(x = reorder(car_name, mean_norm, mean), ymin = mean_norm - se_norm, ymax = mean_norm + se_norm,
                       color = experiment), width = 0.2) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
      scale_y_continuous(name = "CAR Expression (relative to Zeta)") +
      scale_color_manual(values = c("tcsl248" = "#E69F00", "tcsl250" = "#56B4E9"),
                         labels = c("tcsl248" = "TCSL248", "tcsl250" = "TCSL250")) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            legend.position = "right") +
      labs(title = paste0(car_set$title_prefix, " - ", subset, " - CAR Expression (Normalized to Zeta, Mean ± SE)"),
           x = "CAR Construct")
  }

  # Create category-based plots for each subset
  message("Generating category-based plots for CAR set: ", car_set$title_prefix)

  # Generate category-based plots for each subset
  for (subset in cell_subsets) {
    # Differentiation category - non-normalized
    plot_key <- paste0(tolower(subset), "_differentiation_by_category")
    plots[[plot_key]] <- make_category_freq_plot(df_car_set, subset, "Differentiation", normalize = FALSE)

    # Differentiation category - normalized to zeta
    plot_key <- paste0(tolower(subset), "_differentiation_by_category_normalized")
    plots[[plot_key]] <- make_category_freq_plot(df_car_set, subset, "Differentiation", normalize = TRUE)

    # Exhaustion frequency category - non-normalized
    plot_key <- paste0(tolower(subset), "_exhaustion_freq_by_category")
    plots[[plot_key]] <- make_category_freq_plot(df_car_set, subset, "Exhaustion", normalize = FALSE)

    # Exhaustion frequency category - normalized to zeta
    plot_key <- paste0(tolower(subset), "_exhaustion_freq_by_category_normalized")
    plots[[plot_key]] <- make_category_freq_plot(df_car_set, subset, "Exhaustion", normalize = TRUE)

    # Exhaustion MFI category - non-normalized
    plot_key <- paste0(tolower(subset), "_exhaustion_mfi_by_category")
    plots[[plot_key]] <- make_category_mfi_plot(df_car_set, subset, "Exhaustion", normalize = FALSE)

    # Exhaustion MFI category - normalized to zeta
    plot_key <- paste0(tolower(subset), "_exhaustion_mfi_by_category_normalized")
    plots[[plot_key]] <- make_category_mfi_plot(df_car_set, subset, "Exhaustion", normalize = TRUE)

    # Activation category - non-normalized
    plot_key <- paste0(tolower(subset), "_activation_by_category")
    plots[[plot_key]] <- make_category_freq_plot(df_car_set, subset, "Activation", normalize = FALSE)

    # Activation category - normalized to zeta
    plot_key <- paste0(tolower(subset), "_activation_by_category_normalized")
    plots[[plot_key]] <- make_category_freq_plot(df_car_set, subset, "Activation", normalize = TRUE)
  }

  # Save all plots
  message(sprintf("Saving plots for CAR set: %s", car_set$title_prefix))
  for (plot_name in names(plots)) {
    plot_file <- file.path(car_set_dir, paste0(plot_name, ".pdf"))
    ggsave(plot_file, plots[[plot_name]], width = 12, height = 10, create.dir = TRUE)
    message(sprintf("Saving plot: %s", plot_file))
  }
}

message("All processing and plotting complete!")