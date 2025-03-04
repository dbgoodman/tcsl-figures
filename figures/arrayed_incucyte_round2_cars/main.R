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

# Process data
all_data <- process_incucyte_data()

source("figures/arrayed_incucyte_round2_cars/scripts/plot_incucyte.R")

message("Data processing and figure generation complete.") 