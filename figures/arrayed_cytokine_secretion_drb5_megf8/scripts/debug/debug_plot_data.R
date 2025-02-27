library(data.table)
library(here)

# Load data
df_248 <- fread(here('data', 'processed', 'cytokine_secretion', 'tcsl248_processed.csv'))
df_250 <- fread(here('data', 'processed', 'cytokine_secretion', 'tcsl250_processed.csv'))

# Add experiment labels
df_248[, expt := 'TCSL248']
df_250[, expt := 'TCSL250']

# Combine data
df_combined <- rbind(df_248, df_250)

# Define chosen gates
chosen_gate <- c(
  "IL2_pos", "TNFa_pos", "IFNg_pos",
  "IFNg_pos.IL2_pos", "IFNg_pos.TNFa_pos",
  "IL2_pos.TNFa_pos", "IFNg_pos.TNFa_pos.IL2_pos")

# Get plot data
plot_data <- df_combined[
  CAR == T & !(car_name %in% c('POS','NEG')) & !is.na(car_name)][
    gate_clean %in% chosen_gate, .(mean_freq = mean(freq)), 
    by=.(car_name, expt, gate_clean)]

# Print unique car_name values and their levels
cat("\nUnique car_name values:\n")
print(unique(plot_data$car_name))

cat("\nLevels of car_name factor:\n")
print(levels(plot_data$car_name))

# Print first few rows of plot data before and after reordering
cat("\nFirst few rows before reordering:\n")
print(head(plot_data))

plot_data[, scaled_mean_freq := scale(mean_freq), by=.(expt, gate_clean)]
plot_data[, car_name := reorder(car_name, scaled_mean_freq)]

cat("\nFirst few rows after reordering:\n")
print(head(plot_data))

cat("\nLevels of car_name factor after reordering:\n")
print(levels(plot_data$car_name)) 