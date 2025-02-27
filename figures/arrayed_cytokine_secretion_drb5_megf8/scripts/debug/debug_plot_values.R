library(data.table)
library(here)

# Load original data
df_orig <- fread('/Users/dbg/Library/CloudStorage/Box-Box/tcsl/flow/20241112_tcsl248_cytokine/processed_cytokine_data.csv')
df_new <- fread(here('data', 'processed', 'cytokine_secretion', 'tcsl248_processed.csv'))

# Define chosen gates
chosen_gate <- c(
  "IL2_pos", "TNFa_pos", "IFNg_pos",
  "IFNg_pos.IL2_pos", "IFNg_pos.TNFa_pos",
  "IL2_pos.TNFa_pos", "IFNg_pos.TNFa_pos.IL2_pos")

# Function to get plot data
get_plot_data <- function(df, name) {
  plot_data <- df[
    CAR == TRUE & !(car_name %in% c('POS','NEG', '')) & !is.na(car_name)][
      gate_clean %in% chosen_gate, .(mean_freq = mean(freq)), 
      by=.(car_name, subset, gate_clean)]
  
  # Add source label
  plot_data[, source := name]
  return(plot_data)
}

# Get plot data for both datasets
orig_plot <- get_plot_data(df_orig, "original")
new_plot <- get_plot_data(df_new, "new")

# Compare values
cat("\nOriginal plot data for CD3:\n")
print(orig_plot[subset == "CD3"][order(car_name, gate_clean)])

cat("\nNew plot data for CD3:\n")
print(new_plot[subset == "CD3"][order(car_name, gate_clean)])

# Merge and find differences
plot_compare <- merge(
  orig_plot[subset == "CD3"], 
  new_plot[subset == "CD3"], 
  by=c("car_name", "subset", "gate_clean"),
  suffixes=c("_orig", "_new")
)

plot_compare[, diff := mean_freq_orig - mean_freq_new]

cat("\nDifferences in plot values:\n")
print(plot_compare[abs(diff) > 0.01][order(abs(diff), decreasing=TRUE)]) 