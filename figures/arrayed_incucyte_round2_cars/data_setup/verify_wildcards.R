#' Verify that all files matching wildcards in figure_notes.md were copied
#'
#' This script checks Box directories for files matching wildcards and compares
#' with what was actually copied to ensure we didn't miss any files.

library(here)
library(fs)
library(data.table)

# Load configuration
source(here("R", "utils.R"))
config <- load_config()

# Define paths
data_dir <- here("data", "raw", "incucyte_round2")
box_tcsl248_dir <- box_path("incucyte/TCSL248")
box_tcsl250_dir <- box_path("incucyte/TCSL250")

# Function to check if all files matching a pattern were copied
check_pattern_files <- function(box_dir, local_dir, pattern, recursive = FALSE) {
  # List files in Box directory matching the pattern
  box_files <- dir_ls(box_dir, regexp = pattern, recurse = recursive)
  
  # Get just the filenames without the path
  box_filenames <- path_file(box_files)
  
  # List files in local directory
  local_files <- dir_ls(local_dir, recurse = TRUE)
  local_filenames <- path_file(local_files)
  
  # Find files that match the pattern in Box but weren't copied
  missing_files <- box_files[!box_filenames %in% local_filenames]
  
  # Return results
  list(
    pattern = pattern,
    box_dir = box_dir,
    local_dir = local_dir,
    box_files = box_files,
    local_files = local_files,
    missing_files = missing_files
  )
}

# Define patterns to check based on figure_notes.md wildcards
patterns <- list(
  # TCSL248 patterns
  list(
    box_dir = box_tcsl248_dir,
    local_dir = file.path(data_dir, "tcsl248"),
    patterns = c(
      "plate_map.*\\.csv$",
      "car_map.*\\.csv$",
      ".*TCSL248.*\\.txt$",
      ".*TCSL248.*\\.CSV$",
      ".*incucyte.*\\.Rmd$"
    )
  ),
  
  # TCSL250 patterns
  list(
    box_dir = box_tcsl250_dir,
    local_dir = file.path(data_dir, "tcsl250"),
    patterns = c(
      "plate_map.*\\.csv$",
      ".*TCSL250.*\\.txt$",
      ".*TCSL250.*\\.CSV$",
      ".*incucyte.*\\.Rmd$"
    )
  )
)

# Run the checks
all_results <- list()
all_missing_files <- c()

cat("Checking for files matching wildcards in Box that might not have been copied:\n\n")

for (pattern_group in patterns) {
  cat("Checking patterns in", pattern_group$box_dir, ":\n")
  
  for (pattern in pattern_group$patterns) {
    result <- check_pattern_files(pattern_group$box_dir, pattern_group$local_dir, pattern)
    all_results <- c(all_results, list(result))
    
    cat("  Pattern:", pattern, "\n")
    cat("    Files in Box:", length(result$box_files), "\n")
    
    if (length(result$missing_files) > 0) {
      cat("    MISSING FILES:", length(result$missing_files), "\n")
      for (file in result$missing_files) {
        cat("      -", file, "\n")
        all_missing_files <- c(all_missing_files, file)
      }
    } else {
      cat("    All files copied successfully!\n")
    }
    cat("\n")
  }
}

# Summary
cat("\n=== SUMMARY ===\n")
if (length(all_missing_files) > 0) {
  cat("Total missing files:", length(all_missing_files), "\n")
  cat("You may need to update copy_data.R to include these files.\n")
} else {
  cat("All files matching the wildcards have been copied successfully!\n")
  cat("The copy_data.R script appears to be complete.\n")
}

# Check for any additional files in Box that might be relevant
cat("\n=== CHECKING FOR ADDITIONAL POTENTIALLY RELEVANT FILES ===\n")

# Define additional patterns that might be relevant
additional_patterns <- list(
  # TCSL248 additional patterns
  list(
    box_dir = box_tcsl248_dir,
    patterns = c(
      ".*incucyte.*",  # Any file with "incucyte" in the name
      ".*TCSL248.*"    # Any file with "TCSL248" in the name
    )
  ),
  
  # TCSL250 additional patterns
  list(
    box_dir = box_tcsl250_dir,
    patterns = c(
      ".*incucyte.*",  # Any file with "incucyte" in the name
      ".*TCSL250.*"    # Any file with "TCSL250" in the name
    )
  )
)

additional_files <- c()

for (pattern_group in additional_patterns) {
  for (pattern in pattern_group$patterns) {
    box_files <- dir_ls(pattern_group$box_dir, regexp = pattern, recurse = TRUE)
    box_filenames <- path_file(box_files)
    
    # Check if these files were already accounted for
    new_files <- box_files[!box_files %in% unlist(lapply(all_results, function(x) x$box_files))]
    
    if (length(new_files) > 0) {
      additional_files <- c(additional_files, new_files)
    }
  }
}

# Remove duplicates
additional_files <- unique(additional_files)

if (length(additional_files) > 0) {
  cat("Found", length(additional_files), "additional potentially relevant files:\n")
  for (file in additional_files) {
    cat("  +", file, "\n")
  }
  cat("\nConsider whether these files should also be included in copy_data.R\n")
} else {
  cat("No additional potentially relevant files found.\n")
}

# Check if there are any unexpected files in the local directories
cat("\n=== CHECKING FOR UNEXPECTED FILES IN LOCAL DIRECTORIES ===\n")

# Get all files that were expected to be copied
expected_files <- unlist(lapply(all_results, function(x) path_file(x$box_files)))

# Get all files that were actually copied
actual_files <- c(
  path_file(dir_ls(file.path(data_dir, "tcsl248"), recurse = TRUE, type = "file")),
  path_file(dir_ls(file.path(data_dir, "tcsl250"), recurse = TRUE, type = "file"))
)

# Find unexpected files
unexpected_files <- actual_files[!actual_files %in% expected_files]

if (length(unexpected_files) > 0) {
  cat("Found", length(unexpected_files), "unexpected files in local directories:\n")
  for (file in unexpected_files) {
    cat("  ?", file, "\n")
  }
} else {
  cat("No unexpected files found in local directories.\n")
} 