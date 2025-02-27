#' Copy and organize cytokine data from Box structure
#' 
#' This script copies all necessary files for the arrayed cytokine secretion analysis
#' from their original locations in Box to our local data directory.

library(here)
library(fs)

# Load configuration
source(here("R", "utils.R"))
config <- load_config()

# Define source and destination paths
experiments <- list(
  tcsl241_242 = list(
    files = c(
      aurora_combined = "flow/240311_tcsl_242_cytokine/240311_cytokine_both_donors.csv"
    )
  ),
  tcsl248 = list(
    files = c(
      cytokine = "flow/20241112_tcsl248_cytokine/241112_tcsl248_cytokine_TCSL248 Cytokine Unmixed.CSV",
      flowjo_export = "flow/20241112_tcsl248_cytokine/populations.csv",
      plate_map = "flow/20241112_tcsl248_cytokine/plate_map.csv",
      car_map = "flow/20241112_tcsl248_cytokine/car_map.csv"
    )
  ),
  tcsl250 = list(
    files = c(
      cytokine = "flow/20241213.tcsl250 cytokines/2024.12.13 TCSL250 Cytokines_TCSL247 Cytokine Unmixed.CSV",
      flowjo_export = "flow/20241213.tcsl250 cytokines/populations.csv",
      plate_map = "flow/20241213.tcsl250 cytokines/plate_map.csv",
      car_map = "flow/20241213.tcsl250 cytokines/car_map.csv"
    )
  )
)

# Create function to copy files
copy_experiment_files <- function(exp_name, file_list) {
  dest_dir <- here("data", "raw", "cytokine_secretion", exp_name)
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
  for (file_type in names(experiments[[exp_name]]$files)) {
    src_path <- box_path(experiments[[exp_name]]$files[[file_type]])
    dest_path <- path(here("data", "raw", "cytokine_secretion", exp_name), 
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
          here("data", "raw", "cytokine_secretion", "file_manifest.csv"),
          row.names = FALSE) 