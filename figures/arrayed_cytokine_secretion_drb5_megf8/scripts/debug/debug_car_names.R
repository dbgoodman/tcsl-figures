library(data.table)
library(here)

# Load data
df_248 <- fread(here('data', 'processed', 'cytokine_secretion', 'tcsl248_processed.csv'))

# Print unique wells with empty car_names
cat("\nWells with empty car_names:\n")
print(unique(df_248[car_name == "", .(well, car_id)]))

# Print unique car_id values
cat("\nUnique car_id values:\n")
print(unique(df_248$car_id))

# Print unique car_name values
cat("\nUnique car_name values:\n")
print(unique(df_248$car_name))

# Print rows where car_id is empty
cat("\nRows where car_id is empty or NA:\n")
print(df_248[car_id == "" | is.na(car_id), .(well, car_id, car_name, gate_clean, subset)])

# Check unique car_id and car_name mappings
cat("\nUnique car_id to car_name mappings:\n")
print(unique(df_248[, .(car_id, car_name)])) 