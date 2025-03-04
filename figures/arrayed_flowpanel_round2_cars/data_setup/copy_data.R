#' Copy and organize flow cytometry panel data from Box structure
#' 
#' This script copies all necessary files for the flow cytometry panel analysis
#' from their original locations in Box to our local data directory.

library(here)
library(fs)

# Load configuration
source(here("R", "utils.R"))
config <- load_config()

# Define source and destination paths
experiments <- list(
  tcsl248 = list(
    files = c(
      repstim_analysis = "flow/20241116.tcsl248_wk1_repstim/repstim_analysis.Rmd",
      populations = "flow/20241116.tcsl248_wk1_repstim/populations.csv",
      plate_map_1 = "flow/20241116.tcsl248_wk1_repstim/plate_map_1.csv",
      plate_map_2 = "flow/20241116.tcsl248_wk1_repstim/plate_map_2.csv",
      car_map = "flow/20241116.tcsl248_wk1_repstim/car_map.csv"
    )
  ),
  tcsl250 = list(
    files = c(
      populations = "flow/20241219.tcsl250_wk1_repstim/populations.csv",
      plate_map_1 = "flow/20241219.tcsl250_wk1_repstim/plate_map_1.csv",
      plate_map_2 = "flow/20241219.tcsl250_wk1_repstim/plate_map_2.csv",
      car_map = "flow/20241219.tcsl250_wk1_repstim/car_map.csv"
    )
  )
)

# Create function to copy files
copy_experiment_files <- function(exp_name, file_list) {
  dest_dir <- here("data", "raw", "flowpanel_round2", exp_name)
  dir_create(dest_dir)
  
  for (name in names(file_list)) {
    src_path <- box_path(file_list[[name]])
    dest_path <- path(dest_dir, path_file(file_list[[name]]))
    
    # Copy file
    if (file_exists(src_path)) {
      file_copy(src_path, dest_path, overwrite = TRUE)
      message(sprintf("Copied %s to %s", src_path, dest_path))
    } else {
      warning(sprintf("Source file not found: %s", src_path))
    }
  }
}

# Copy all files
for (exp_name in names(experiments)) {
  copy_experiment_files(exp_name, experiments[[exp_name]]$files)
}

# Create a manifest of copied files
manifest <- data.frame(
  experiment = character(),
  file_type = character(),
  source_path = character(),
  dest_path = character(),
  stringsAsFactors = FALSE
)

for (exp_name in names(experiments)) {
  # Add files to manifest
  for (file_type in names(experiments[[exp_name]]$files)) {
    src_path <- box_path(experiments[[exp_name]]$files[[file_type]])
    dest_path <- path(here("data", "raw", "flowpanel_round2", exp_name), 
                     path_file(experiments[[exp_name]]$files[[file_type]]))
    
    manifest <- rbind(manifest, data.frame(
      experiment = exp_name,
      file_type = file_type,
      source_path = src_path,
      dest_path = dest_path,
      stringsAsFactors = FALSE
    ))
  }
}

# Save manifest
write.csv(manifest, 
          here("data", "raw", "flowpanel_round2", "file_manifest.csv"),
          row.names = FALSE) 