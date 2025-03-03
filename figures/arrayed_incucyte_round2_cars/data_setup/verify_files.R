#' Verify that all expected files were copied correctly
#'
#' This script checks that all the wildcards in the copy_data.R script
#' match the actual files that were copied.

library(here)
library(fs)
library(data.table)

# Define paths
data_dir <- here("data", "raw", "incucyte_round2")

# Check TCSL248 files
tcsl248_dir <- file.path(data_dir, "tcsl248")
tcsl248_raw_dir <- file.path(tcsl248_dir, "raw_data")

# Expected files for TCSL248
expected_tcsl248_files <- c(
  # Metadata files
  file.path(tcsl248_dir, "car_map.csv"),
  file.path(tcsl248_dir, "plate_map_wk1.csv"),
  file.path(tcsl248_dir, "plate_map_wk2.csv"),
  file.path(tcsl248_dir, "incucyte_analysis.Rmd"),
  
  # Raw data files
  file.path(tcsl248_raw_dir, "20241118_TCSL248_wk0_incucyte_count.txt"),
  file.path(tcsl248_raw_dir, "20241118_TCSL248_wk0_incucyte_total_area.txt"),
  file.path(tcsl248_raw_dir, "20241118_TCSL248_wk0_incucyte_total_int_intensity.txt"),
  file.path(tcsl248_raw_dir, "20241121_TCSL248_wk1_incucyte_count.txt"),
  file.path(tcsl248_raw_dir, "20241121_TCSL248_wk1_incucyte_total_area.txt"),
  file.path(tcsl248_raw_dir, "20241121_TCSL248_wk1_incucyte_total_int_intensity.txt"),
  file.path(tcsl248_raw_dir, "20241115_TCSL248_wk1count_redo_Count Raw plate2.CSV")
)

# Check TCSL250 files
tcsl250_dir <- file.path(data_dir, "tcsl250")
tcsl250_raw_dir <- file.path(tcsl250_dir, "raw_data")

# Expected files for TCSL250
expected_tcsl250_files <- c(
  # Metadata files
  file.path(tcsl250_dir, "plate_map_r2.csv"),
  file.path(tcsl250_dir, "plate_map_r8.csv"),
  file.path(tcsl250_dir, "plate_map_r24.csv"),
  file.path(tcsl250_dir, "plate_map_purified-r28.csv"),
  file.path(tcsl250_dir, "incucyte_analysis.Rmd"),
  
  # Raw data files
  file.path(tcsl250_raw_dir, "20241118_TCSL250_wk0_incucyte_count_r8.txt"),
  file.path(tcsl250_raw_dir, "20241118_TCSL250_wk0_incucyte_total_area_r8.txt"),
  file.path(tcsl250_raw_dir, "20241118_TCSL250_wk0_incucyte_total_int_intensity_r8.txt"),
  file.path(tcsl250_raw_dir, "20241118_TCSL250_wk0_incucyte_count_r24.txt"),
  file.path(tcsl250_raw_dir, "20241118_TCSL250_wk0_incucyte_total_area_r24.txt"),
  file.path(tcsl250_raw_dir, "20241118_TCSL250_wk0_incucyte_total_int_intensity_r24.txt"),
  file.path(tcsl250_raw_dir, "20241126_TCSL250_wk1_incucyte_count_r2.txt"),
  file.path(tcsl250_raw_dir, "20241126_TCSL250_wk1_incucyte_total_area_r2.txt"),
  file.path(tcsl250_raw_dir, "20241126_TCSL250_wk1_incucyte_total_int_intensity_r2.txt"),
  file.path(tcsl250_raw_dir, "20241223_TCSL250_wk0_incucyte_count_purified-r28.txt"),
  file.path(tcsl250_raw_dir, "20241223_TCSL250_wk0_incucyte_total_area_purified-r28.txt"),
  file.path(tcsl250_raw_dir, "20241223_TCSL250_wk0_incucyte_total_int_intensity_purified-r28.txt")
)

# Combine all expected files
all_expected_files <- c(expected_tcsl248_files, expected_tcsl250_files)

# Check if each expected file exists
missing_files <- c()
for (file in all_expected_files) {
  if (!file_exists(file)) {
    missing_files <- c(missing_files, file)
  }
}

# Print results
cat("Verification of copied files:\n")
cat("Total expected files:", length(all_expected_files), "\n")

if (length(missing_files) == 0) {
  cat("All expected files were copied successfully!\n")
} else {
  cat("Missing files:", length(missing_files), "\n")
  for (file in missing_files) {
    cat(" -", file, "\n")
  }
}

# Check if there are any unexpected files
all_actual_files <- c(
  dir_ls(tcsl248_dir, recurse = TRUE, type = "file"),
  dir_ls(tcsl250_dir, recurse = TRUE, type = "file")
)

unexpected_files <- setdiff(all_actual_files, all_expected_files)

if (length(unexpected_files) > 0) {
  cat("\nUnexpected files found:", length(unexpected_files), "\n")
  for (file in unexpected_files) {
    cat(" +", file, "\n")
  }
} else {
  cat("\nNo unexpected files found.\n")
}

# Check if the data directories exist and have the right structure
cat("\nDirectory structure verification:\n")
dirs_to_check <- c(
  data_dir,
  tcsl248_dir,
  tcsl248_raw_dir,
  tcsl250_dir,
  tcsl250_raw_dir
)

for (dir in dirs_to_check) {
  if (dir_exists(dir)) {
    cat(" ✓", dir, "exists\n")
  } else {
    cat(" ✗", dir, "does not exist\n")
  }
}

# Check the content of a few key files to ensure they're not empty
cat("\nChecking content of key files:\n")
key_files <- c(
  file.path(tcsl248_dir, "car_map.csv"),
  file.path(tcsl248_dir, "plate_map_wk1.csv"),
  file.path(tcsl250_dir, "plate_map_r2.csv")
)

for (file in key_files) {
  if (file_exists(file)) {
    file_size <- file_size(file)
    if (file_size > 0) {
      cat(" ✓", file, "has content (", file_size, "bytes)\n")
      
      # Try to read the first few lines
      tryCatch({
        content <- readLines(file, n = 3)
        cat("   First few lines:\n")
        for (line in content) {
          cat("   ", line, "\n")
        }
      }, error = function(e) {
        cat("   Error reading file:", e$message, "\n")
      })
    } else {
      cat(" ✗", file, "is empty\n")
    }
  } else {
    cat(" ✗", file, "does not exist\n")
  }
} 