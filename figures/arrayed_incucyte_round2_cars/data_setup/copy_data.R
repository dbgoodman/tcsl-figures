#' Copy and organize incucyte data from Box structure
#' 
#' This script copies all necessary files for the arrayed incucyte killing analysis
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
      car_map = "incucyte/TCSL248/car_map.csv",
      plate_map_wk1 = "incucyte/TCSL248/plate_map_wk1.csv",
      plate_map_wk2 = "incucyte/TCSL248/plate_map_wk2.csv",
      incucyte_analysis = "incucyte/TCSL248/incucyte_analysis.Rmd"
    ),
    raw_data_files = c(
      # Week 0 data
      count_wk0 = "incucyte/TCSL248/20241118_TCSL248_wk0_incucyte_count.txt",
      area_wk0 = "incucyte/TCSL248/20241118_TCSL248_wk0_incucyte_total_area.txt",
      intensity_wk0 = "incucyte/TCSL248/20241118_TCSL248_wk0_incucyte_total_int_intensity.txt",
      # Week 1 data
      count_wk1 = "incucyte/TCSL248/20241121_TCSL248_wk1_incucyte_count.txt",
      area_wk1 = "incucyte/TCSL248/20241121_TCSL248_wk1_incucyte_total_area.txt",
      intensity_wk1 = "incucyte/TCSL248/20241121_TCSL248_wk1_incucyte_total_int_intensity.txt",
      # Additional data
      count_wk1_redo = "incucyte/TCSL248/20241115_TCSL248_wk1count_redo_Count Raw plate2.CSV"
    )
  ),
  tcsl250 = list(
    files = c(
      # Plate maps for different replicates
      plate_map_r2 = "incucyte/TCSL250/plate_map_r2.csv",
      plate_map_r8 = "incucyte/TCSL250/plate_map_r8.csv",
      plate_map_r24 = "incucyte/TCSL250/plate_map_r24.csv",
      plate_map_purified_r28 = "incucyte/TCSL250/plate_map_purified-r28.csv",
      incucyte_analysis = "incucyte/TCSL250/incucyte_analysis.Rmd"
    ),
    raw_data_files = c(
      # Week 0 data - no replicate specified (missing files)
      count_wk0 = "incucyte/TCSL250/20241118_TCSL250_wk0_incucyte_count.txt",
      area_wk0 = "incucyte/TCSL250/20241118_TCSL250_wk0_incucyte_total_area.txt",
      
      # Week 0 data - r8 replicate
      count_wk0_r8 = "incucyte/TCSL250/20241118_TCSL250_wk0_incucyte_count_r8.txt",
      area_wk0_r8 = "incucyte/TCSL250/20241118_TCSL250_wk0_incucyte_total_area_r8.txt",
      intensity_wk0_r8 = "incucyte/TCSL250/20241118_TCSL250_wk0_incucyte_total_int_intensity_r8.txt",
      
      # Week 0 data - r24 replicate
      count_wk0_r24 = "incucyte/TCSL250/20241118_TCSL250_wk0_incucyte_count_r24.txt",
      area_wk0_r24 = "incucyte/TCSL250/20241118_TCSL250_wk0_incucyte_total_area_r24.txt",
      intensity_wk0_r24 = "incucyte/TCSL250/20241118_TCSL250_wk0_incucyte_total_int_intensity_r24.txt",
      
      # Week 1 data - r2 replicate
      count_wk1_r2 = "incucyte/TCSL250/20241126_TCSL250_wk1_incucyte_count_r2.txt",
      area_wk1_r2 = "incucyte/TCSL250/20241126_TCSL250_wk1_incucyte_total_area_r2.txt",
      intensity_wk1_r2 = "incucyte/TCSL250/20241126_TCSL250_wk1_incucyte_total_int_intensity_r2.txt",
      
      # Purified data - r28 replicate
      count_purified_r28 = "incucyte/TCSL250/20241223_TCSL250_wk0_incucyte_count_purified-r28.txt",
      area_purified_r28 = "incucyte/TCSL250/20241223_TCSL250_wk0_incucyte_total_area_purified-r28.txt",
      intensity_purified_r28 = "incucyte/TCSL250/20241223_TCSL250_wk0_incucyte_total_int_intensity_purified-r28.txt"
    )
  )
)

# Create function to copy files
copy_experiment_files <- function(exp_name, file_list) {
  dest_dir <- here("data", "raw", "incucyte_round2", exp_name)
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

# Copy raw data files
copy_raw_data_files <- function(exp_name, file_list) {
  dest_dir <- here("data", "raw", "incucyte_round2", exp_name, "raw_data")
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
  copy_raw_data_files(exp_name, experiments[[exp_name]]$raw_data_files)
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
  # Add regular files to manifest
  for (file_type in names(experiments[[exp_name]]$files)) {
    src_path <- box_path(experiments[[exp_name]]$files[[file_type]])
    dest_path <- path(here("data", "raw", "incucyte_round2", exp_name), 
                     path_file(experiments[[exp_name]]$files[[file_type]]))
    
    manifest <- rbind(manifest, data.frame(
      experiment = exp_name,
      file_type = file_type,
      source_path = src_path,
      dest_path = dest_path,
      stringsAsFactors = FALSE
    ))
  }
  
  # Add raw data files to manifest
  for (file_type in names(experiments[[exp_name]]$raw_data_files)) {
    src_path <- box_path(experiments[[exp_name]]$raw_data_files[[file_type]])
    dest_path <- path(here("data", "raw", "incucyte_round2", exp_name, "raw_data"), 
                     path_file(experiments[[exp_name]]$raw_data_files[[file_type]]))
    
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
          here("data", "raw", "incucyte_round2", "file_manifest.csv"),
          row.names = FALSE) 