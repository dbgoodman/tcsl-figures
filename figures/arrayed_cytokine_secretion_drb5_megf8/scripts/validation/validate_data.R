library(data.table)
library(here)

# Load both datasets
df_orig <- fread('/Users/dbg/Library/CloudStorage/Box-Box/tcsl/flow/20241112_tcsl248_cytokine/processed_cytokine_data.csv')
df_new <- fread(here('data', 'processed', 'cytokine_secretion', 'tcsl248_processed.csv'))

# Basic comparison
cat("\nDimensions:")
cat("\nOriginal:", dim(df_orig))
cat("\nNew:", dim(df_new))

# Compare column names
col_diff <- setdiff(names(df_orig), names(df_new))
if (length(col_diff) > 0) {
  cat("\n\nColumns in original but not in new:", col_diff)
}
col_diff <- setdiff(names(df_new), names(df_orig))
if (length(col_diff) > 0) {
  cat("\nColumns in new but not in original:", col_diff)
}

# Compare values for CD3 subset
cat("\n\nComparing CD3 values:")
orig_cd3 <- df_orig[subset == "CD3" & CAR == TRUE & !(car_name %in% c("POS", "NEG", "")), 
                    .(mean_freq = mean(freq)), by=.(car_name, gate_clean)]
new_cd3 <- df_new[subset == "CD3" & CAR == TRUE & !(car_name %in% c("POS", "NEG", "")), 
                  .(mean_freq = mean(freq)), by=.(car_name, gate_clean)]

# Merge and compare
cd3_compare <- merge(orig_cd3, new_cd3, by=c("car_name", "gate_clean"), suffixes=c("_orig", "_new"))
cd3_compare[, diff := mean_freq_orig - mean_freq_new]

cat("\nSummary of differences:")
print(cd3_compare[abs(diff) > 0.01])

# Compare unique values in key columns
cat("\n\nUnique values comparison:")
for (col in c("car_name", "gate_clean", "subset")) {
  cat("\n", col, ":")
  cat("\nOriginal:", paste(sort(unique(df_orig[[col]])), collapse=", "))
  cat("\nNew:", paste(sort(unique(df_new[[col]])), collapse=", "))
}

# Compare raw values for a specific example
cat("\n\nDetailed comparison for DRB5 CD3 IL2_pos:")
print(df_orig[car_name == "DRB5" & subset == "CD3" & gate_clean == "IL2_pos", 
             .(car_name, subset, gate_clean, freq)])
print(df_new[car_name == "DRB5" & subset == "CD3" & gate_clean == "IL2_pos", 
             .(car_name, subset, gate_clean, freq)]) 