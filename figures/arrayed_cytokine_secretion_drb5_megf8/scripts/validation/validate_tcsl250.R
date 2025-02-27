library(data.table)
library(here)

# Load our processed data
df_new <- fread(here('data', 'processed', 'cytokine_secretion', 'tcsl250_processed.csv'))

# Load reference processed data from Box
df_ref <- fread('/Users/dbg/Library/CloudStorage/Box-Box/tcsl/flow/20241213.tcsl250 cytokines/processed_cytokine_data.csv')

# Basic validation
cat("\nDimensions check:\n")
cat("Reference:", dim(df_ref), "\n")
cat("New:", dim(df_new), "\n")

# Column comparison
cat("\nColumns in reference but not in new:\n")
print(setdiff(names(df_ref), names(df_new)))
cat("\nColumns in new but not in reference:\n")
print(setdiff(names(df_new), names(df_ref)))

# Compare key statistics for CD8+ CAR+ cells
get_summary <- function(df) {
  df[CAR == TRUE & subset == "CD8" & 
     gate_clean %in% c("IL2_pos", "IFNg_pos", "TNFa_pos"),
     .(mean_freq = mean(freq),
       sd_freq = sd(freq)), 
     by = .(gate_clean, car_name)]
}

cat("\nReference summary:\n")
print(get_summary(df_ref))
cat("\nNew summary:\n")
print(get_summary(df_new))

# Check exact matches for key columns
check_exact_match <- function(col) {
  ref_vals <- sort(unique(df_ref[[col]]))
  new_vals <- sort(unique(df_new[[col]]))
  cat("\nUnique", col, "values match:", identical(ref_vals, new_vals))
  if (!identical(ref_vals, new_vals)) {
    cat("\nIn reference but not in new:", setdiff(ref_vals, new_vals))
    cat("\nIn new but not in reference:", setdiff(new_vals, ref_vals))
  }
}

check_exact_match("gate_clean")
check_exact_match("car_name")
check_exact_match("subset")

# Compare frequencies for key cytokines
cat("\nComparing frequencies for key cytokines (CD8+ CAR+ cells):\n")
key_gates <- c("IL2_pos", "IFNg_pos", "TNFa_pos")
for (gate in key_gates) {
  cat("\nGate:", gate, "\n")
  ref_freqs <- df_ref[CAR == TRUE & subset == "CD8" & gate_clean == gate, 
                     .(car_name, freq)]
  new_freqs <- df_new[CAR == TRUE & subset == "CD8" & gate_clean == gate, 
                     .(car_name, freq)]
  
  # Sort both by car_name and freq to ensure matching order
  ref_freqs <- ref_freqs[order(car_name, freq)]
  new_freqs <- new_freqs[order(car_name, freq)]
  
  # Check if frequencies match exactly
  freq_match <- all.equal(ref_freqs$freq, new_freqs$freq)
  cat("Frequencies match exactly:", identical(freq_match, TRUE), "\n")
  if (!identical(freq_match, TRUE)) {
    cat("Differences found:", freq_match, "\n")
  }
} 