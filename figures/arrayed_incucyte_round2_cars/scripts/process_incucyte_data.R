#' Process incucyte data for TCSL248 and TCSL250 experiments
#'
#' This script loads and processes the incucyte data for the arrayed CAR experiments,
#' following the approach used in the original analysis Rmd files.

library(data.table)
library(ggplot2)
library(patchwork)
library(stringr)
library(here)
library(fs)

# Define paths
data_dir <- here("data", "raw", "incucyte_round2")
processed_dir <- here("data", "processed", "incucyte_round2")
dir_create(processed_dir)

# Function to load and process a plate map (from reference files)
load_plate_metadata <- function(plate_map_filepath, exclude_empty = TRUE, rows = LETTERS[1:8], columns = 12) {
  # Load plate map
  plate_map <- fread(plate_map_filepath, header = FALSE)
  
  # Truncate or pad rows
  if (nrow(plate_map) > length(rows)) {
    plate_map <- plate_map[1:length(rows), ]
  } else if (nrow(plate_map) < length(rows)) {
    plate_map <- rbind(plate_map, matrix(NA, nrow = length(rows) - nrow(plate_map), ncol = ncol(plate_map)))
  }
  
  # Truncate or pad columns
  if (ncol(plate_map) > columns) {
    plate_map <- plate_map[, 1:columns, with = FALSE]
  } else if (ncol(plate_map) < columns) {
    plate_map <- cbind(plate_map, matrix(NA, nrow = nrow(plate_map), ncol = columns - ncol(plate_map)))
  }
  
  # Add row labels
  plate_map[, rn := rows[seq_len(nrow(plate_map))]]  # Add row names as a column
  
  # Set column names
  setnames(plate_map, c(as.character(seq_len(columns)), "rn"))
  
  # Convert to long format
  long_data <- melt(as.data.table(plate_map), id.vars = "rn", variable.name = "col", value.name = "name")
  long_data[, `:=`(well = paste0(rn, col), row = rn, col = as.integer(col))]
  
  # Exclude empty wells if required
  if (exclude_empty) {
    long_data <- long_data[!is.na(name)]
  }
  
  # Return the long format data
  return(long_data[, .(well, row, col, name)])
}

# Function to load a single measurement file
load.single.measurement <- function(data.path) {
  # Ensure we're using absolute paths
  if (!file.exists(data.path)) {
    stop(paste("File not found:", data.path))
  }
  
  data.dt <- fread(data.path)
  data.dt <- data.dt[, 2:ncol(data.dt)]
  data.melt.dt <- melt(data.dt, id.vars='Elapsed', variable.name='test', value.name='value')
  data.melt.dt <- data.melt.dt[, c('well','image') := tstrsplit(test, ", Image")]
  data.melt.dt[, value := as.numeric(value)]
  return(data.melt.dt)
}

# Process TCSL248 data
process_tcsl248 <- function() {
  message("Processing TCSL248 data...")
  
  # Define paths
  tcsl248_dir <- file.path(data_dir, "tcsl248")
  tcsl248_raw_dir <- file.path(tcsl248_dir, "raw_data")
  
  # List data files
  tcsl248.filepaths.dt <- data.table(
    data.path = list.files(path = tcsl248_raw_dir, pattern = '.+TCSL248.+\\.txt$', full.names = TRUE)
  )
  
  tcsl248.filepaths.dt[,
    `:=`(
      measurement = gsub('.+wk\\d_(\\w+)\\.txt', '\\1', data.path),
      week = gsub('.+wk(\\d)_\\w+\\.txt', '\\1', data.path)
    )]
  
  # Load plate maps
  platemap.dt <- rbind(
    load_plate_metadata(file.path(tcsl248_dir, "plate_map_wk1.csv"))[, week := 0],
    load_plate_metadata(file.path(tcsl248_dir, "plate_map_wk2.csv"))[, week := 1]
  )
  
  # Split `name` column into `car_id` and `ratio` at the `.`
  platemap.dt[, c("car_id", "ratio") := tstrsplit(name, ".", fixed=TRUE)]
  platemap.dt[car_id %in% c('U','K'), ratio := NA_character_]

  # Load and combine data
  combined.data.dt <- tcsl248.filepaths.dt[,
    load.single.measurement(data.path), 
    by=c('measurement','week')]
  
  # Merge with plate map
  combined.data.dt <- merge(combined.data.dt[, week := as.numeric(week)], platemap.dt, by=c('well','week'))
  
  # Load car map
  car_map <- fread(file.path(tcsl248_dir, "car_map.csv"))
  car_map[, car_id := factor(car_id)]
  car_map[, car_name := factor(car_name)]
  
  # Merge with car map
  combined.data.dt <- car_map[combined.data.dt, on='car_id']
  combined.data.dt[is.na(car_name), car_name := car_id]
  combined.data.dt[, car_id := factor(car_id, levels=c(1:16,'U','K'))]
  
  # Filter out extreme outlier values
  # First, calculate reasonable bounds for each well/measurement/week
  combined.data.dt[, `:=`(
    median_value = median(value, na.rm = TRUE),
    iqr_value = IQR(value, na.rm = TRUE)
  ), by = .(well, measurement, week)]
  
  # Define outlier threshold - using a much more lenient threshold (10 times IQR instead of 5)
  combined.data.dt[, `:=`(
    lower_bound = median_value - 30 * iqr_value,
    upper_bound = median_value + 30 * iqr_value
  )]

  # Flag outliers, but NEVER flag Elapsed == 0 values
  combined.data.dt[, is_outlier := value < lower_bound | value > upper_bound]
  combined.data.dt[Elapsed == 0, is_outlier := FALSE]  # Protect Elapsed == 0 values

  # Replace outliers with NA
  combined.data.dt[is_outlier == TRUE, value := NA]
  
  # Remove temporary columns used for outlier detection
  combined.data.dt[, c("median_value", "iqr_value", "lower_bound", "upper_bound", "is_outlier") := NULL]
  
  # Normalize data
  combined.data.dt[, value_norm := value/(value[Elapsed == 0]+0.001), 
    by=c('well','measurement', 'week')]
  combined.data.dt[, mean_value_norm := mean(value_norm, na.rm=TRUE),
    by=c('ratio','car_name','Elapsed','measurement','week')]
  combined.data.dt[, std_dev := sd(value_norm, na.rm=TRUE),
    by=c('ratio','car_name','Elapsed','measurement','week')]
  combined.data.dt[, std_error := std_dev/sqrt(.N),
    by=c('ratio','car_name','Elapsed','measurement','week')]
  
  # Normalize to control
  combined.data.dt[week == 0, ctrl_val := mean(mean_value_norm[car_name == 'U'], na.rm=TRUE), by=c('Elapsed','measurement')]
  combined.data.dt[week == 1, ctrl_val := mean(mean_value_norm[car_name == 'K'], na.rm=TRUE), by=c('Elapsed','measurement')]
  
  # Add normalization to no T cells
  combined.data.dt[, mean_value_none := mean_value_norm/ctrl_val,
    by=c('Elapsed','measurement','week')]
  combined.data.dt[, value_none := value_norm/ctrl_val,
    by=c('ratio','Elapsed','measurement','week')]
  combined.data.dt[, std_dev_value_none := sd(value_none, na.rm=TRUE),
    by=c('ratio','car_name','Elapsed','measurement','week')]
  combined.data.dt[, std_err_value_none := std_dev_value_none/sqrt(.N),
    by=c('ratio','car_name','Elapsed','measurement','week')]  
  # Add experiment identifier
  combined.data.dt[, experiment := "tcsl248"]
  
  return(combined.data.dt)
}

# Process TCSL250 data
process_tcsl250 <- function() {
  message("Processing TCSL250 data...")
  
  # Define paths
  tcsl250_dir <- file.path(data_dir, "tcsl250")
  tcsl250_raw_dir <- file.path(tcsl250_dir, "raw_data")
  
  # List data files
  filepaths.dt <- data.table(
    data.path = list.files(path = tcsl250_raw_dir, pattern = '.+TCSL250.+\\.txt$', full.names = TRUE)
  )
  
  # Extract measurement type, week, and plate from filenames
  filepaths.dt[, `:=`(
    measurement = gsub('.+wk\\d_(\\w+)_(purified-)?r\\d+\\.txt', '\\1', data.path),
    week = gsub('.+wk(\\d)_\\w+_(purified-)?r\\d+\\.txt', '\\1', data.path),
    plate = gsub('.+wk\\d_\\w+_((purified-)?r\\d+)\\.txt', '\\1', data.path)
  )]
  
  # Load plate maps
  platemap.dt <- rbind(
    load_plate_metadata(file.path(tcsl250_dir, "plate_map_r8.csv"))[, `:=`(week = 0, plate = 'r8')],
    load_plate_metadata(file.path(tcsl250_dir, "plate_map_r2.csv"))[, `:=`(week = 1, plate = 'r2')],
    load_plate_metadata(file.path(tcsl250_dir, "plate_map_purified-r28.csv"))[, `:=`(week = 0, plate = 'purified-r28')],
    load_plate_metadata(file.path(tcsl250_dir, "plate_map_r24.csv"))[, `:=`(week = 0, plate = 'r24')]
  )
  
  # Split name column into car_id and ratio
  platemap.dt[, c("car_id", "ratio") := tstrsplit(name, ".", fixed=TRUE)]
  platemap.dt[car_id %in% c('U','K'), ratio := NA_character_]
  # Load and combine data
  combined.data.dt <- filepaths.dt[,
    load.single.measurement(data.path), 
    by=c('measurement', 'week', 'plate')]
  
  # Merge with plate map
  combined.data.dt <- merge(
    combined.data.dt[, week := as.numeric(week)], 
    platemap.dt, 
    by=c('well', 'week', 'plate')
  )

  
  # Load car map (using TCSL248's car map as mentioned in the notes)
  car_map <- fread(file.path(data_dir, "tcsl248", "car_map.csv"))
  car_map[, car_id := factor(car_id)]
  car_map[, car_name := factor(car_name)]
  
  # Merge with car map
  combined.data.dt <- car_map[combined.data.dt, on='car_id']
  combined.data.dt[is.na(car_name), car_name := car_id]
  combined.data.dt[, car_id := factor(car_id, levels=c(1:16,'U','K'))]
  
  
  # Add experiment identifier
  combined.data.dt[, experiment := "tcsl250"]
  
  # make purified it's own week
  combined.data.dt[, week := as.factor(week)]
  combined.data.dt[plate == 'purified-r28', week := paste0(week, "_purified")]
  
  # print(combined.data.dt[, .N, by=.(week, plate, ratio, measurement)])
  # Remove temporary columns
  combined.data.dt[, c("name") := NULL]
  combined.data.dt[, purified := grepl("purified", week)]

  #  remove col A from purified plate
  combined.data.dt <- combined.data.dt[plate != 'purified-r28' | col != 1]

  # combine individual images for purified plate 
  combined.data.dt <- rbind(
    combined.data.dt[plate != 'purified-r28'],
    combined.data.dt[plate == 'purified-r28', 
      list(value=mean(value, na.rm=TRUE), image=well[1], test=well[1]), 
      by=setdiff(names(combined.data.dt), c('value', 'image', 'test'))]
  )
  
  # Filter out extreme outlier values
  # First, calculate reasonable bounds for each well/measurement/week/plate
  combined.data.dt[, `:=`(
    median_value = median(value, na.rm = TRUE),
    iqr_value = IQR(value, na.rm = TRUE)
  ), by = .(well, measurement, week, plate)]
  
  # Define outlier threshold - using a much more lenient threshold (10 times IQR instead of 5)
  combined.data.dt[, `:=`(
    lower_bound = median_value - 30 * iqr_value,
    upper_bound = median_value + 30 * iqr_value
  )]
  
  # Flag outliers, but NEVER flag Elapsed == 0 values
  combined.data.dt[, is_outlier := value < lower_bound | value > upper_bound]
  combined.data.dt[Elapsed == 0, is_outlier := FALSE]  # Protect Elapsed == 0 values
  
  # Replace outliers with NA
  combined.data.dt[is_outlier == TRUE, value := NA]
  
  # Remove temporary columns used for outlier detection
  combined.data.dt[, c("median_value", "iqr_value", "lower_bound", "upper_bound", "is_outlier") := NULL]
  
  # Normalize data
  combined.data.dt[, value_norm := value/(value[Elapsed == 0]+0.001), 
    by=c('well','measurement', 'week', 'plate')]
  combined.data.dt[, mean_value_norm := mean(value_norm, na.rm=TRUE),
    by=c('ratio','car_name','Elapsed','measurement','week')]
  combined.data.dt[, std_dev := sd(value_norm, na.rm=TRUE),
    by=c('ratio','car_name','Elapsed','measurement','week')]
  combined.data.dt[, std_error := std_dev/sqrt(.N),
    by=c('ratio','car_name','Elapsed','measurement','week')]
  
  # Normalize to control
  combined.data.dt[, ctrl_val := mean(mean_value_norm[car_name == 'U'], na.rm=TRUE), 
                  by=c('Elapsed','measurement','week')]
  
  # Add normalization to no T cells
  combined.data.dt[, mean_value_none := mean_value_norm/ctrl_val,
    by=c('Elapsed','measurement','week')]
  combined.data.dt[, value_none := value_norm/ctrl_val,
    by=c('ratio','Elapsed','measurement','week')]
  combined.data.dt[, std_dev_value_none := sd(value_none, na.rm=TRUE),
    by=c('ratio','car_name','Elapsed','measurement','week')]
  combined.data.dt[, std_err_value_none := std_dev_value_none/sqrt(.N),
    by=c('ratio','car_name','Elapsed','measurement','week')]
  
  # Apply reasonable bounds to value_none and value_norm
  combined.data.dt[value_none < -2, value_none := NA]
  combined.data.dt[value_none > 5, value_none := NA]
  combined.data.dt[value_norm < -2, value_norm := NA]
  combined.data.dt[value_norm > 5, value_norm := NA]
  
  return(combined.data.dt)
}

# Process all incucyte data
process_incucyte_data <- function() {
  message("Processing all incucyte data...")
  
  # Process TCSL248 data
  tcsl248_data <- process_tcsl248()
  
  # Process TCSL250 data
  tcsl250_data <- process_tcsl250()
  
  # Combine all data
  all_data <- rbindlist(list(tcsl248_data, tcsl250_data), use.names = TRUE, fill = TRUE)
  
  message("Data processing complete.")
  
  # Return the processed data
  return(all_data)
}

# Execute the processing
all_data <- process_incucyte_data()

# Print summary
cat("Processed incucyte data summary:\n")
cat("Number of experiments:", length(unique(all_data$experiment)), "\n")
cat("Number of replicates:", length(unique(all_data$replicate)), "\n")
cat("Number of time points:", length(unique(all_data$Elapsed)), "\n")
cat("Number of wells:", length(unique(all_data$well)), "\n")
cat("Number of CARs:", length(unique(all_data$car_id)), "\n")

# Return the processed data
all_data

# Save the processed data to CSV
fwrite(all_data, file.path(processed_dir, "all_incucyte_data.csv")) 