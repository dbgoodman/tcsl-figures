#' Amino Acid Sequence Alignment Demo
#' 
#' This script demonstrates the use of the amino acid sequence alignment functions
#' to visualize alignments between different CAR variants.

library(ggplot2)
library(data.table)
library(grid)
library(gridExtra)

# Source required files
source("figures/shared/utils.R")
source("figures/shared/car_sets.R")
source("figures/shared/car_colors.R")
source("figures/aa_seq_alignment/scripts/alignment_functions.R")

# Create output directory if it doesn't exist
figure_dir <- "figures/aa_seq_alignment/output"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Load CAR data
car_data <- load_car_data()

# Print available car sets for debugging
cat("Available car sets:\n")
for (set_name in names(car_sets)) {
  cat(sprintf("- %s: %s\n", set_name, paste(car_sets[[set_name]]$cars, collapse=", ")))
}

# Set plot dimensions - wider to accommodate labels
plot_width <- 8
plot_height <- 6

# Define the alignments to generate
alignments <- list(
  list(name = "megf8", ref_car = "MEGF8", car_set = car_sets[["megf8_all_variants"]]),
  list(name = "cd28", ref_car = "CD28", car_set = car_sets[["cd28_variants"]]),
  list(name = "tnr9", ref_car = "TNR9", car_set = car_sets[["tnr9_variants"]]),
  list(name = "drb5", ref_car = "DRB5", car_set = car_sets[["drb5_all_variants"]])
)

# Generate and save each alignment plot
for (alignment in alignments) {
  # Generate the plot
  plot <- aa_seq_alignment_plot(alignment$car_set, alignment$ref_car, car_data)
  
  # Save as PDF
  filename <- file.path(figure_dir, paste0(alignment$name, "_alignment.pdf"))
  ggsave(filename, plot, width = plot_width, height = plot_height)
}

# Print completion message
cat("Amino acid sequence alignment plots have been saved to", figure_dir, "\n") 