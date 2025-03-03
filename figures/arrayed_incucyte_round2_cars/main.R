#' Main script for generating figures for arrayed incucyte round 2 CAR experiments
#'
#' This script sources the necessary scripts to process the data and generate figures.

# Load required libraries
library(data.table)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(viridis)

# Create figure directory if it doesn't exist
figure_dir <- "figures/arrayed_incucyte_round2_cars/output"
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

# Clean up any existing output files
unlink(figure_dir, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

# Source data processing and plotting scripts
source("figures/arrayed_incucyte_round2_cars/scripts/process_incucyte_data.R")
source("figures/arrayed_incucyte_round2_cars/scripts/plot_incucyte.R")

# Process data
all_data <- process_incucyte_data()

# Define final timepoints map
final_timepoints_map <- list(
  tcsl248 = list(
    week_0 = 60,  # Week 0: 60 hours
    week_1 = 80   # Week 1: 80 hours
  ),
  tcsl250 = list(
    week_0 = 64,
    week_1 = 64
  )
)



# Plot TCSL248 time course
plot_tcsl248_time_course(all_data, figure_dir)

# Plot TCSL250 time course
plot_tcsl250_time_course(all_data, figure_dir)

# Plot comparison
plot_comparison(all_data, figure_dir)

# Plot heatmap
plot_heatmap(all_data, figure_dir)

message("Data processing and figure generation complete.") 