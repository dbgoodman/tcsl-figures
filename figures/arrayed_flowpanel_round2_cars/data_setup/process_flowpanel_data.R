#' Process flow cytometry panel data
#' 
#' This script processes flow cytometry panel data from TCSL248 and TCSL250 experiments.
#' It loads, cleans, and prepares the data for plotting.

library(data.table)
library(ggplot2)
library(patchwork)
library(here)

#' Process flow cytometry data for a specific experiment
#'
#' @param experiment_id The experiment ID (e.g., "tcsl248" or "tcsl250")
#' @param data_dir The directory containing the raw data files
#' @param rep_metadata A data.table with replicate metadata (columns: rep, rep_type, stimmed, rep_lbl)
#' @return A processed data.table with all the necessary annotations
process_flow_data <- function(experiment_id, data_dir = NULL, rep_metadata = NULL) {
  
  # Set data directory if not provided
  if (is.null(data_dir)) {
    data_dir <- here("data", "raw", "flowpanel_round2", experiment_id)
  }
  
  # Load raw data
  df_raw <- fread(file.path(data_dir, "populations.csv"))
  
  # Keep *.fcs only
  df_raw <- df_raw[grepl('.fcs', V1)]
  
  # Extract plate and well
  df_raw[, `:=`(plate = as.integer(sub(".*Plate_([0-9]+).*", "\\1", V1)),
            well = sub(".*-([A-Z][0-9]+) .*", "\\1", V1))]
  
  # Remove unneeded columns
  df_raw[, grep("^V", names(df_raw)) := NULL]
  
  # Melt
  df_melt <- unique(melt(df_raw, id.vars=c('plate', 'well'), variable.name='gate_str', value.name='value'))
  
  # Remove freq parent from live car
  df_melt <- df_melt[!(grepl('Live \\w+ \\| Freq. of Parent', gate_str))]
  
  # Extract columns
  df_melt[, `:=`(
    gate_str = sub(' \\|.*', '', gate_str),
    gate = sub(" \\|.*", "", sub(".*/", "", gate_str)),  # Remove everything after '|' and handle top-level gates
    parent_gate = fifelse(grepl("^([^/]+) \\|", gate_str),        # Top-level gates (no preceding '/')
                     NA_character_,                         # Set parent to NA for top-level gates
                     sub(".*/([^/]+)/[^/]+ \\|.*", "\\1", gate_str)),  # Extract parent for nested gates
    measure = gsub('.*\\| ', '', gate_str)  # else extract measure from gate_str and process
  )]
  
  df_melt[grepl('Freq', measure), measure := 'freq']
  df_melt[grepl('Count', measure), measure := 'count']
  df_melt[, measure := gsub('Geometric Mean \\((.*)\\)', 'gmean_\\1', measure)]
  df_melt[, measure := gsub('Mean \\((.*)\\)', 'mean_\\1', measure)]
  
  df_melt[grepl('mean', measure), parent_gate := gate]
  df_melt[grepl('mean', measure), gate := gsub('.*mean_(.*)', '\\1',  measure)]
  df_melt[grepl('mean', measure), gate_str := paste0(gate_str, '/', gsub('.*mean_(.*)', '\\1',  measure))]
  df_melt[grepl('mean', measure), measure := gsub('(.*mean)_.*', '\\1_mfi',  measure)]
  
  # If parent is 'Single Cells', make parent 'Lymphocytes' and keep only gate str with 'Freq. of Lymphocytes'
  df_melt[parent_gate == 'Single Cells', parent_gate := 'Lymphocytes']
  df_melt <- df_melt[is.na(parent_gate) | parent_gate != 'Lymphocytes' | (parent_gate == 'Lymphocytes' & measure == 'freq')]
  
  # Add CAR status
  df_melt[, CAR := grepl('CAR', gate_str)]
  
  # Add CD4/CD8 status
  df_melt[, subset := 'CD3']
  df_melt[grepl('/(CD[48])/', gate_str), subset := gsub('.*/(CD[48])/.*', '\\1', gate_str)]
  
  # Get parent subset and car status
  df_melt[, parent_subset := subset]
  df_melt[grepl('^CD[48]$', parent_gate), parent_subset := 'CD3']
  df_melt[, parent_CAR := CAR]
  df_melt[gate == 'Live CAR', parent_CAR := FALSE]
  
  # Separate freq and count
  df_melt <- dcast(df_melt, ... ~ measure, value.var='value')
  
  # Remove count/freq from gate_str
  df_melt[, gate_str := sub(' \\|.*', '', gate_str)]
  
  # Consolidate separate freq/count rows
  df_melt[is.na(parent_gate), freq := freq[!is.na(freq)], by= .(plate, well, gate, CAR, subset)]
  df_melt[is.na(parent_gate), count := count[!is.na(count)], by= .(plate, well, gate, CAR, subset)]
  df_melt <- unique(df_melt)
  
  # Set NA freq to 0
  df_melt[is.na(freq), freq := 0]
  df_melt[, count := as.numeric(count)]
  df_melt[, freq := as.numeric(freq)]
  df_melt[, gmean_mfi := as.numeric(gmean_mfi)]
  df_melt[, mean_mfi := as.numeric(mean_mfi)]
  
  # Propagate counts iteratively
  while (any(is.na(df_melt[!(gate %in% c('beads','Lymphocytes')), count]))) {
    # Get new counts
    df_new_counts <- df_melt[!is.na(count), .(CAR, gate, subset, well, plate, count)][df_melt[is.na(count)],   
      on=.(CAR=parent_CAR, gate=parent_gate, subset=parent_subset, well, plate)][
        !is.na(count)][, 
        .(plate, well, gate=i.gate, CAR=i.CAR, subset=i.subset,
          count=round(as.numeric(freq)/100 * as.numeric(count)))]
    
    # Add new counts onto df_melt
    df_melt[df_new_counts, 
            on = .(plate, well, gate, CAR, subset), 
            count := fifelse(is.na(count), i.count, count)]
    
  }
  
  # Clean gates
  # Create a copy of df_melt to work on gate names
  df_melt[, gate_clean := gate]
  
  # Apply regex transformations
  df_melt[, gate_clean := gsub("\\+\\s*", "pos", gate_clean)]    # Replace '+' with 'pos' and handle trailing spaces
  df_melt[, gate_clean := gsub("\\-\\s*", "neg", gate_clean)]    # Replace '-' with 'neg' and handle trailing spaces
  df_melt[, gate_clean := gsub(",\\s*", "", gate_clean)]         # Remove spaces and commas after ','
  df_melt[, gate_clean := gsub("\\s+", "", gate_clean)]          # Remove all remaining spaces
  
  # Simplify Q gates
  df_melt[, gate_clean := gsub("^Q[0-9]+:CD", "CD", gate_clean)]  # Remove 'QX:' prefix and keep CDs intact for clarity
  
  # Case adjustments
  df_melt[, gate_clean := gsub("Lymphocytes", "lymphocytes", gate_clean)]  # Lowercase lymphocytes
  df_melt[, gate_clean := gsub("LiveCAR", "live_car", gate_clean)]         # Replace Live CAR with live_car
  df_melt[, gate_clean := gsub("LiveUntr", "live_untr", gate_clean)]       # Replace Live Untr with live_untr
  
  # Optional: Convert to factor and reorder levels if needed
  df_melt[, gate_clean := factor(gate_clean, levels = unique(gate_clean))]
  
  # Apply gsub to rename each memory subset
  df_melt[, gate_clean := gsub("CD45ROnegCD62Lpos", "Naive", gate_clean)]
  df_melt[, gate_clean := gsub("CD45ROposCD62Lpos", "Mem", gate_clean)]
  df_melt[, gate_clean := gsub("CD45ROposCD62Lneg", "Eff", gate_clean)]
  df_melt[, gate_clean := gsub("CD45ROnegCD62Lneg", "Emra", gate_clean)]
  
  # Convert to a factor to preserve ordering
  df_melt[, gate_clean := factor(gate_clean, levels = unique(gate_clean))]
  
  # Get count per bead
  df_melt[, count_per_bead := count / count[gate == 'beads'][1], by=.(plate, well)]
  
  # Combine subgates for 25, 27 and 127
  # Identify the gates to combine
  combine_gates <- c("CD25", "CD27", "CD127")
  
  # Loop through each gate and combine its subgates
  for (gate in combine_gates) {
    # Find all subgates for the current gate
    subgate_pattern <- paste0(gate, 'pos')  # Match gates starting with the gate name
    subgates <- grep(subgate_pattern, df_melt$gate_clean, value = TRUE)
    
    # Check if any subgates are present (skip if not found)
    if (length(subgates) > 0) {
      # Sum counts for all subgates and assign to the combined gate
      combined_counts <- df_melt[gate_clean %in% subgates, .(
        count = sum(count, na.rm = TRUE), 
        freq = sum(freq, na.rm= TRUE),
        parent_gate, parent_CAR, parent_subset
      ), 
      by = .(plate, well, CAR, subset)]
      
      combined_counts[, gate_clean := gate]  # Assign the combined gate name
      
      # Add the combined gate to the main table (keep original gates)
      df_melt <- rbind(df_melt, combined_counts, fill = TRUE)
    }
  }
  
  # Convert gate_clean to factor again
  df_melt[, gate_clean := factor(gate_clean, levels = unique(gate_clean))]
  
  # Create proper CD3 data by aggregating CD4 and CD8 data
  # First, get all unique gate_clean values that have MFI data
  mfi_gates <- unique(df_melt[!is.na(mean_mfi), gate_clean])
  
  # For each gate_clean and plate/well/CAR combination, calculate CD3 values
  cd3_data <- data.table()
  
  for (gc in unique(df_melt$gate_clean)) {
    # Get CD4 and CD8 data for this gate
    cd4_cd8_data <- df_melt[subset %in% c("CD4", "CD8") & gate_clean == gc]
    
    if (nrow(cd4_cd8_data) > 0) {
      # Calculate total counts and weighted MFI values
      cd3_agg <- cd4_cd8_data[, .(
        count = sum(count, na.rm = TRUE),
        # Calculate weighted mean for MFI values (weighted by count)
        mean_mfi = weighted.mean(mean_mfi, count, na.rm = TRUE),
        gmean_mfi = weighted.mean(gmean_mfi, count, na.rm = TRUE),
        # Calculate weighted frequency (weighted by count)
        freq = weighted.mean(freq, count, na.rm = TRUE),
        # Store parent information
        parent_gate = "CD3",
        parent_subset = "CD3",
        parent_CAR = unique(parent_CAR)[1],
        # Store other needed columns
        gate_str = paste0("CD3/", gc),
        gate = gc
      ), by = .(plate, well, CAR, gate_clean)]
      
      # Set subset to CD3
      cd3_agg[, subset := "CD3"]
      
      # Add to CD3 data
      cd3_data <- rbind(cd3_data, cd3_agg, fill = TRUE)
    }
  }
  
  # Add the CD3 data to the main data
  df_melt <- rbind(df_melt, cd3_data, fill = TRUE)
  
  # Add plate map
  # Function to load and process a plate map
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
  
  # Load plate maps
  plate_map_1 <- load_plate_metadata(file.path(data_dir, "plate_map_1.csv"))
  plate_map_2 <- load_plate_metadata(file.path(data_dir, "plate_map_2.csv"))
  
  # Add plate identifiers to the plate maps
  plate_map_1[, plate := 1]
  plate_map_2[, plate := 2]
  
  # Combine plate maps
  plate_map <- rbind(plate_map_1, plate_map_2)
  
  # Map plate map names onto df_melt
  df_melt <- merge(df_melt, plate_map, by = c("plate", "well"), all.x = TRUE)
  
  # Rename the `name` field to `well_id`
  setnames(df_melt, "name", "well_id")
  
  # Split `well_id` into `rep` and `car_id`
  df_melt[, `:=`(
    rep = substr(well_id, 1, 1),  # Extract the first character as `rep`
    car_id = substr(well_id, 2, nchar(well_id))  # Extract the remainder as `car_id`
  )]
  
  # Convert `car_id` to a factor with correct order (1-16, U)
  df_melt[, car_id := factor(car_id, levels = c(as.character(1:16), "U"))]
  
  # Join `df_melt` with `car_map.csv` on `car_id`
  car_map <- fread(file.path(data_dir, "car_map.csv"))  # Load the mapping file
  car_map[, car_id := factor(car_id)]
  car_map[, car_name := factor(car_name)]
  df_melt <- merge(df_melt, car_map, by = "car_id", all.x = TRUE)
  df_melt[is.na(car_name), car_name := car_id]  # Replace missing car names with car_id
  
  # Fix any issues with the rep column
  if(!is.character(df_melt$rep)) {
    # Create a safe character version based on plate
    df_melt[, rep := as.character(plate)]
  }
  
  # Apply replicate metadata if provided
  if (!is.null(rep_metadata)) {
    # Merge replicate metadata with df_melt
    df_melt <- merge(df_melt, rep_metadata, by = "rep", all.x = TRUE)
  } else {
    # Default replicate metadata (for backward compatibility)
    df_melt[rep == 'X', rep_type := 'Unstim']
    df_melt[rep %in% c('A','B'), rep_type := '1:8']
    df_melt[rep == 'C', rep_type := '1:4']
    df_melt[, stimmed := rep_type != 'Unstim']
    df_melt[, rep_lbl := factor(rep, levels=c('X','A','B','C'), labels = c('Unstim', '1:8a', '1:8b', '1:4'))]
  }
  
  # Assign gate categories
  df_melt[gate_clean %in% c('CD25','CD27','CD127'), gate_category := 'Activation']
  df_melt[gate_clean %in% c('Eff','Emra','Mem','Naive'), gate_category := 'Differentiation']
  df_melt[gate_clean %in% c('TIM3','LAG3','CD39','PD1'), gate_category := 'Exhaustion']
  
  return(df_melt)
}

#' Function to process data for a specific experiment and save it
#'
#' @param experiment_id The experiment ID (e.g., "tcsl248" or "tcsl250")
#' @param rep_metadata A data.table with replicate metadata (columns: rep, rep_type, stimmed, rep_lbl)
#' @param output_dir Directory to save the processed data
#' @return A processed data.table
process_and_save_flow_data <- function(experiment_id, rep_metadata = NULL, output_dir = NULL) {
  # Process the data
  df_processed <- process_flow_data(experiment_id, rep_metadata = rep_metadata)
  
  # Set output directory if not provided
  if (is.null(output_dir)) {
    output_dir <- here("data", "processed", "flowpanel_round2")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Save the processed data
  output_file <- file.path(output_dir, paste0(experiment_id, "_processed.rds"))
  saveRDS(df_processed, output_file)
  
  message(sprintf("Processed data saved to %s", output_file))
  
  return(df_processed)
} 