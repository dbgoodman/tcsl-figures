#' Process a single flow cytometry experiment
#'
#' This script provides a function to process a single flow cytometry experiment
#' (TCSL248 or TCSL250) and generate plots.

library(here)
library(data.table)

#' Process a single flow cytometry experiment
#'
#' @param experiment_id The experiment ID (e.g., "tcsl248" or "tcsl250")
#' @param rep_metadata A data.table with replicate metadata
#' @param output_dir Directory to save the processed data
#' @param plots_dir Directory to save the plots
#' @param bead_count_file Path to the bead count file (optional)
#' @return The processed data.table
process_experiment <- function(experiment_id, rep_metadata, 
                              output_dir = NULL, plots_dir = NULL, 
                              bead_count_file = NULL) {
  # Set default output directory if not provided
  if (is.null(output_dir)) {
    output_dir <- here("data", "processed", "flowpanel_round2")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Set default plots directory if not provided
  if (is.null(plots_dir)) {
    plots_dir <- here("figures", "arrayed_flowpanel_round2_cars", "plots", experiment_id)
    dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Process the data
  message(sprintf("Processing %s flow cytometry data...", experiment_id))
  df_processed <- process_flow_data(experiment_id, rep_metadata = rep_metadata)
  
  # Save the processed data as CSV instead of RDS
  output_file <- file.path(output_dir, paste0(experiment_id, "_processed.csv"))
  fwrite(df_processed, output_file)
  message(sprintf("Processed data saved to %s", output_file))
  
  # Generate plots
  message(sprintf("Generating plots for %s...", experiment_id))
  save_all_plots(df_processed, experiment_id, bead_count_file, plots_dir)
  message(sprintf("Plots saved to %s", plots_dir))
  
  # Return the processed data
  return(df_processed)
} 