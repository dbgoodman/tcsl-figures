# Test script for cytokine data processing

library(here)
library(data.table)

# Source the processing function
source(here("figures", "arrayed_cytokine_secretion_drb5_megf8", "process_cytokine_data.R"))

# Process TCSL248 data
tcsl248_data <- process_experiment_data(
  experiment_dir = here("data", "raw", "cytokine_secretion", "tcsl248"),
  output_file = here("data", "processed", "cytokine_secretion", "tcsl248_processed.csv")
)

# Basic validation checks
cat("\nBasic validation of processed data:\n")
cat("Number of rows:", nrow(tcsl248_data), "\n")
cat("Number of unique gates:", length(unique(tcsl248_data$gate_clean)), "\n")
cat("\nUnique gates:\n")
print(sort(unique(tcsl248_data$gate_clean)))
cat("\nUnique CAR names:\n")
print(sort(unique(tcsl248_data$car_name)))
cat("\nSample of processed data:\n")
print(tcsl248_data[CAR == TRUE & gate_clean %in% c("IL2_pos", "IFNg_pos") & 
                   subset == "CD8" & !is.na(car_name)][1:5])

# Compare with some key metrics from original analysis
cat("\nSummary of cytokine positive cells for CD8+ CAR+ cells:\n")
print(tcsl248_data[
  CAR == TRUE & 
  subset == "CD8" & 
  gate_clean %in% c("IL2_pos", "IFNg_pos", "TNFa_pos") &
  !is.na(car_name) & 
  !(car_name %in% c("POS", "NEG")),
  .(mean_freq = mean(freq), sd_freq = sd(freq)), 
  by = .(gate_clean)]) 