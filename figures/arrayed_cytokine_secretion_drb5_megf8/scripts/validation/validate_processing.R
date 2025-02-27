# Validate that our processed data matches the original

library(data.table)
library(here)

# Load both datasets
original <- fread('/Users/dbg/Library/CloudStorage/Box-Box/tcsl/flow/20241112_tcsl248_cytokine/processed_cytokine_data.csv')
new <- fread(here('data', 'processed', 'cytokine_secretion', 'tcsl248_processed.csv'))

# Basic comparisons
cat("\nBasic dataset comparisons:\n")
cat("Original dimensions:", dim(original), "\n")
cat("New dimensions:", dim(new), "\n")

# Compare column names
cat("\nColumn differences:\n")
print(setdiff(names(original), names(new)))
print(setdiff(names(new), names(original)))

# Compare key statistics for cytokine frequencies
cat("\nComparing cytokine frequencies for CD8+ CAR+ cells:\n")

get_summary <- function(data) {
  data[CAR == TRUE & 
       subset == "CD8" & 
       gate_clean %in% c("IL2_pos", "IFNg_pos", "TNFa_pos") &
       !is.na(car_name) & 
       !(car_name %in% c("POS", "NEG")),
       .(mean_freq = mean(freq), sd_freq = sd(freq)), 
       by = .(gate_clean, car_name)]
}

orig_summary <- get_summary(original)
new_summary <- get_summary(new)

# Compare summaries
cat("\nOriginal summary:\n")
print(orig_summary)
cat("\nNew summary:\n")
print(new_summary)

# Check for exact matches in key columns
check_exact_match <- function(col) {
  orig_sorted <- sort(unique(original[[col]]))
  new_sorted <- sort(unique(new[[col]]))
  cat("\nComparing unique values in", col, ":\n")
  cat("Original:", length(orig_sorted), "values\n")
  cat("New:", length(new_sorted), "values\n")
  if (!identical(orig_sorted, new_sorted)) {
    cat("Differences found:\n")
    cat("In original but not new:", setdiff(orig_sorted, new_sorted), "\n")
    cat("In new but not original:", setdiff(new_sorted, orig_sorted), "\n")
  } else {
    cat("Exact match!\n")
  }
}

# Check key columns
key_cols <- c("gate_clean", "car_name", "subset")
for (col in key_cols) {
  check_exact_match(col)
} 