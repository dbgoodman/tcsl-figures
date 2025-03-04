#' Copy reference Rmarkdown files for adaptation
#' 
library(here)
library(fs)

# Load configuration
source(here("R", "utils.R"))
config <- load_config()

# Define source files
reference_files <- c(
  tcsl248 = "flow/20241116.tcsl248_wk1_repstim/repstim_analysis.Rmd"
  # Note: TCSL250 doesn't have a repstim_analysis.Rmd file according to the figure_notes.md
)

# Copy files with meaningful names
dest_dir <- here("figures", "arrayed_flowpanel_round2_cars", "reference")
dir_create(dest_dir)

for (exp in names(reference_files)) {
  src_path <- box_path(reference_files[[exp]])
  dest_path <- path(dest_dir, sprintf("%s_analysis.Rmd", exp))
  
  if (file_exists(src_path)) {
    file_copy(src_path, dest_path, overwrite = TRUE)
    message(sprintf("Copied %s to %s", src_path, dest_path))
  } else {
    warning(sprintf("Source file not found: %s", src_path))
  }
} 