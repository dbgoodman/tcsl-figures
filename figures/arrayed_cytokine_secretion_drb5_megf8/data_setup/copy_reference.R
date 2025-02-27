#' Copy reference Rmarkdown files for adaptation
#' 
library(here)
library(fs)

# Load configuration
source(here("R", "utils.R"))
config <- load_config()

# Define source files
reference_files <- c(
  tcsl248 = "flow/20241112_tcsl248_cytokine/cytokine_analysis.Rmd",
  tcsl250 = "flow/20241213.tcsl250 cytokines/cytokine_analysis.Rmd"
)

# Copy files with meaningful names
dest_dir <- here("figures", "arrayed_cytokine_secretion_drb5_megf8", "reference")
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