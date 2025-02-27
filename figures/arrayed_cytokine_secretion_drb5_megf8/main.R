# Master script for cytokine secretion analysis
library(here)

# Process the raw data
source(here('figures', 'arrayed_cytokine_secretion_drb5_megf8', 'process_cytokine_data.R'))

# Generate all plots
source(here('figures', 'arrayed_cytokine_secretion_drb5_megf8', 'plot_cytokines.R'))

message("All processing and plotting complete!") 