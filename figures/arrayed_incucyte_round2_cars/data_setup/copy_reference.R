#' Copy reference Rmarkdown files for adaptation
#' 
library(here)
library(fs)

# Load configuration
source(here("R", "utils.R"))
config <- load_config()

# Define source files
reference_files <- c(
  tcsl248 = "incucyte/TCSL248/incucyte_analysis.Rmd",
  tcsl250 = "incucyte/TCSL250/incucyte_analysis.Rmd"
)

# Copy files with meaningful names
dest_dir <- here("figures", "arrayed_incucyte_round2_cars", "reference")
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