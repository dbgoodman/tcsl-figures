#' Process cytokine data from flow cytometry experiments
#' 
#' This script contains functions to process cytokine data from flow cytometry experiments,
#' based on the analysis from TCSL248 and TCSL250.

library(data.table)
library(stringr)
library(here)
library(dplyr)

#' Process flow cytometry population data for a single experiment
#' 
#' @param experiment_dir Directory containing the experiment files
#' @param output_file Path to save the processed data
#' @return Processed data.table with all cytokine populations
process_experiment_data <- function(experiment_dir, output_file = NULL) {
  
  # Load raw population data
  df_raw <- fread(file.path(experiment_dir, 'populations.csv'))
  
  # keep *.fcs only
  df_raw <- df_raw[grepl('.fcs', V1)]
  
  # Extract plate and well
  df_raw[, `:=`(plate = as.integer(sub(".*Plate_([0-9]+).*", "\\1", V1)),
            well = sub(".*-([A-Z][0-9]+) .*", "\\1", V1))]
  
  # remove unneeded columns
  df_raw[, grep("^V", names(df_raw)) := NULL]
  
  # melt
  df_melt <- unique(melt(df_raw, id.vars=c('plate', 'well'), 
                        variable.name='gate_str', value.name='value'))
  
  # extract columns
  df_melt[, `:=`(
    gate = sub(" \\|.*", "", sub(".*/", "", gate_str)),  
    parent_gate = fifelse(grepl("^([^/]+) \\|", gate_str),        
                   NA_character_,                         
                   gsub(".*/([^/]+)/[^/]+ \\|.*", "\\1", 
                        gsub('Q10: IFNg\\+ , TNF-a\\+/', '', gate_str))),  
    measure = ifelse(grepl("\\| Freq", gate_str), "freq", 
                    ifelse(grepl("\\| Count", gate_str), "count", NA))
  )]
  
  # if parent is 'Single Cells', make parent 'Lymphocytes'
  df_melt[parent_gate == 'Single Cells', parent_gate := 'Lymphocytes']
  df_melt <- df_melt[is.na(parent_gate) | 
                     parent_gate != 'Lymphocytes' | 
                     grepl('Freq. of Lymphocytes', gate_str)]
  
  # add CAR status
  df_melt[, CAR := grepl('_car', gate_str)]
  
  # add CD4/CD8 status
  df_melt[, subset := 'None']
  df_melt[grepl('/(CD[348])/', gate_str), 
          subset := gsub('.*/(CD[348])/.*', '\\1', gate_str)]
  
  # get parent subset and car status
  df_melt[, parent_subset := subset]
  df_melt[grepl('^CD[348]$', parent_gate), parent_subset := 'None']
  df_melt[, parent_CAR := CAR]
  df_melt[gate == 'live_car', parent_CAR := FALSE]
  
  # separate freq and count
  df_melt <- dcast(df_melt, ... ~ measure, value.var='value')
  
  # remove count/freq from gate_str
  df_melt[, gate_str := sub(' \\|.*', '', gate_str)]
  
  # consolidate separate freq/count rows
  df_melt[is.na(parent_gate), 
          freq := freq[!is.na(freq)], 
          by= .(plate, well, gate, CAR, subset)]
  df_melt[is.na(parent_gate), 
          count := count[!is.na(count)], 
          by= .(plate, well, gate, CAR, subset)]
  df_melt <- unique(df_melt)
  
  # set NA freq to 0
  df_melt[is.na(freq), freq := 0]
  
  # Propagate counts iteratively
  while (any(is.na(df_melt[!(gate %in% c('beads','Lymphocytes')), count]))) {
    df_new_counts <- df_melt[!is.na(count), 
                            .(CAR, gate, subset, well, plate, count)][
      df_melt[is.na(count)],   
      on=.(CAR=parent_CAR, 
          gate=parent_gate, 
          subset=parent_subset, 
          well, plate)][
            !is.na(count)][, 
                          .(plate, well, gate=i.gate, CAR=i.CAR, subset=i.subset,
                            count=round(freq/100 * count))]
    
    df_melt[df_new_counts, 
            on = .(plate, well, gate, CAR, subset), 
            count := fifelse(is.na(count), i.count, count)]
  }
  
  # Clean gate names
  df_melt[, gate_clean := gate]
  df_melt[, gate_clean := gsub("TNF-a", "TNFa", gate_clean)]
  df_melt[, gate_clean := gsub("\\+", "_pos", gate_clean)]
  df_melt[, gate_clean := gsub("\\-", "_neg", gate_clean)]
  df_melt[, gate_clean := gsub(" , ", ".", gate_clean)]
  df_melt[, gate_clean := gsub("^\\s+", "", gate_clean)]
  df_melt[, gate_clean := gsub("^Q[0-9]+: ", "", gate_clean)]
  
  # Handle triple cytokine gate exactly as in reference
  df_melt[grepl('Q10: IFNg\\+ , TNF-a\\+/', gate_str), 
    gate_clean := paste0('IFNg_pos.TNFa_pos.', gate_clean)]
  
  # Convert to factor to preserve ordering
  df_melt[, gate_clean := factor(gate_clean, levels = unique(gate_clean))]
  
  # Load plate map and car map exactly as in reference
  plate_map <- fread(file.path(experiment_dir, "plate_map.csv"), header = FALSE, na.strings = "")
  
  # Truncate or pad rows to 8
  if (nrow(plate_map) > 8) {
    plate_map <- plate_map[1:8, ]
  } else if (nrow(plate_map) < 8) {
    plate_map <- rbind(plate_map, matrix(NA, nrow = 8 - nrow(plate_map), ncol = ncol(plate_map)))
  }
  
  # Truncate or pad columns to 12
  if (ncol(plate_map) > 12) {
    plate_map <- plate_map[, 1:12, with = FALSE]
  } else if (ncol(plate_map) < 12) {
    plate_map <- cbind(plate_map, matrix(NA, nrow = nrow(plate_map), ncol = 12 - ncol(plate_map)))
  }
  
  # Add row labels and set column names
  plate_map[, rn := LETTERS[1:8]]
  setnames(plate_map, c(as.character(1:12), "rn"))
  
  # Convert to long format
  long_plate <- melt(plate_map, id.vars = "rn", variable.name = "col", value.name = "name")
  long_plate[, `:=`(
    well = paste0(rn, col),
    row = rn,
    col = as.integer(col)
  )]
  
  # Exclude empty wells
  long_plate <- long_plate[!is.na(name)]
  
  # Merge plate map with data - keep all rows from df_melt
  df_melt <- merge(df_melt, long_plate[, .(well, row, col, name)], by = "well", all.x = TRUE)
  
  # Rename name to car_id
  setnames(df_melt, "name", "car_id")
  
  # Convert car_id to factor with correct order (1-16, NEG, POS)
  df_melt[, car_id := factor(car_id, levels = c(as.character(1:16), "NEG", "POS"))]
  
  # Load and merge CAR map
  car_map <- fread(file.path(experiment_dir, "car_map.csv"))
  car_map[, car_id := factor(car_id)]
  car_map[, car_name := factor(car_name)]
  df_melt <- merge(df_melt, car_map, by = "car_id", all.x = TRUE)
  df_melt[is.na(car_name), car_name := as.character(car_id)]  # Replace missing car names with car_id
  
  # Save if output file specified
  if (!is.null(output_file)) {
    fwrite(df_melt, output_file)
  }
  
  return(df_melt)
}

#' Process TCSL241/242 data which has a simpler format
#' 
#' @param input_file Path to the combined donor data CSV
#' @param output_file Path to save the processed data
#' @return Processed data.table in format compatible with other experiments
process_241_242_data <- function(input_file, output_file = NULL) {
  # Load data
  dt <- fread(input_file)
  
  # Create mapping for consistent CAR names, using padded numbers
  car_name_map <- c(
    '01' = 'CD28',
    '02' = 'TNR9',
    '05' = 'Zeta',
    '06' = 'MEGF8',
    '12' = 'DRB5'
  )
  
  # Process data
  dt <- melt(dt, id.vars=c('donor', 'sample'))[
    # Pad single digit numbers with 0
    , sample := sprintf("%02d", as.numeric(sample))][
    # Clean variable names and make consistent with 248/250
    , variable := gsub('_pct', '', variable)][
    , variable := dplyr::case_when(
        variable == "ifng" ~ "IFNg",
        variable == "tnfa" ~ "TNFa",
        variable == "il2" ~ "IL2",
        TRUE ~ variable
    )][
    # Convert to same format as 248/250
    , gate_clean := paste0(variable, '_pos')
  ]
  
  # Fix experiment names to be uppercase
  dt[, expt := toupper(donor)]
  
  # Add CAR column and subset
  dt[, `:=`(
    CAR = TRUE,
    subset = 'CD3',
    car_name = sample
  )]
  
  # Map CAR names to be consistent
  dt[car_name %in% names(car_name_map), 
     car_name := car_name_map[car_name]]
  
  # Remove unwanted samples and those not in the mapping
  dt <- dt[car_name %in% car_name_map]
  
  # Rename columns to match 248/250 format
  setnames(dt, 'value', 'freq')
  
  # Save if output file specified
  if (!is.null(output_file)) {
    fwrite(dt, output_file)
  }
  
  return(dt)
}

#' Process TCSL247 data which has a different format and multiple donors
#' 
#' @param input_dir Directory containing the TCSL247 data files
#' @param output_file Path to save the processed data
#' @return Processed data.table in format compatible with other experiments
process_247_data <- function(input_dir, output_file = NULL) {
  # Load data
  dt <- fread(file.path(input_dir, "TCSL247 Cytokine Unmixed.csv"))
  plate_map <- fread(file.path(input_dir, "plate_map.csv"), header=FALSE)
  sample_metadata <- fread(file.path(input_dir, "sample_metadata.csv"))
  
  # Extract relevant columns (% Parent for cytokines and cell count)
  cytokine_cols <- c(
    "IL2" = "IL2 % Parent",
    "IFNg" = "IFNg % Parent",
    "TNFa" = "TNFa % Parent"
  )
  
  # Select and rename columns
  dt_subset <- dt[, c("V1", cytokine_cols, "cells Count"), with=FALSE]
  setnames(dt_subset, 
           old=c("V1", cytokine_cols, "cells Count"),
           new=c("well", names(cytokine_cols), "cell_count"))
  
  # Clean well IDs - extract just the well position (e.g., A1, B2)
  dt_subset[, well := sub(".*-([A-Z][0-9]+) Well.*", "\\1", well)]
  
  # Melt to long format
  dt_long <- melt(dt_subset, 
                  id.vars=c("well", "cell_count"),
                  variable.name="cytokine",
                  value.name="freq")
  
  # Add gate_clean column in same format as other experiments
  dt_long[, gate_clean := paste0(cytokine, "_pos")]
  
  # Process sample metadata
  setnames(sample_metadata, "test/invivo", "test")
  sample_metadata[, `:=`(
    donor = case_when(
      test == "V" ~ "D1V",
      donor == 1 ~ "D1",
      donor == 2 ~ "D2"
    ),
    car_name = car
  )]
  
  # Map DRB5.CD28 to DRB528C
  sample_metadata[car_name == "DRB5.CD28", car_name := "DRB528C"]
  
  # Create mapping from well to sample/donor
  # First create a data.table with all well positions
  wells <- CJ(row=LETTERS[1:nrow(plate_map)], col=1:ncol(plate_map))
  wells[, well := paste0(row, col)]
  
  # Add sample numbers from plate map
  wells[, sample := as.vector(t(as.matrix(plate_map)))]
  
  # Debug: Print wells before filtering
  print("Wells before filtering:")
  print(wells[well == "A11"])
  
  # Keep only numeric samples and convert to integer, but handle single-digit numbers
  wells <- wells[!is.na(sample) & sample != ""]
  wells[, sample := as.integer(gsub("^0+", "", sample))]
  
  # Debug: Print wells after filtering
  print("Wells after filtering:")
  print(wells[well == "A11"])
  
  # Merge well map with sample metadata
  well_metadata <- merge(wells[, .(well, sample)], 
                        sample_metadata[, .(number, donor, car_name)],
                        by.x="sample", by.y="number", all.x=TRUE)
  
  # Debug: Print well metadata
  print("Well metadata after merge:")
  print(well_metadata[well == "A11"])
  
  # Merge metadata with data
  dt_final <- merge(dt_long, well_metadata, by="well")
  
  # Debug: Print final data
  print("Final data for A11:")
  print(dt_final[well == "A11"])
  
  # Add experiment name for each donor
  dt_final[, expt := paste0("TCSL247_", donor)]
  
  # Add required columns to match other experiments
  dt_final[, `:=`(
    CAR = car_name != "KO",  # Set CAR to FALSE for knockout controls
    subset = "CD3"
  )]
  
  # Save if output file specified
  if (!is.null(output_file)) {
    fwrite(dt_final, output_file)
  }
  
  return(dt_final)
}

# Process all experiments
if (interactive() || !exists('TESTING')) {
  # Process TCSL248 and TCSL250
  df_248 <- process_experiment_data(
    here('data', 'raw', 'cytokine_secretion', 'tcsl248'),
    here('data', 'processed', 'cytokine_secretion', 'tcsl248_processed.csv')
  )
  
  df_250 <- process_experiment_data(
    here('data', 'raw', 'cytokine_secretion', 'tcsl250'),
    here('data', 'processed', 'cytokine_secretion', 'tcsl250_processed.csv')
  )
  
  # Process TCSL241/242
  df_241_242 <- process_241_242_data(
    here('data', 'raw', 'cytokine_secretion', 'tcsl241_242', '240311_cytokine_both_donors.csv'),
    here('data', 'processed', 'cytokine_secretion', 'tcsl241_242_processed.csv')
  )
  
  # Process TCSL247
  df_247 <- process_247_data(
    here('data', 'raw', 'cytokine_secretion', 'tcsl247'),
    here('data', 'processed', 'cytokine_secretion', 'tcsl247_processed.csv')
  )
}

# Example usage:
# process_experiment_data("data/raw/cytokine_secretion/tcsl248",
#                        "data/processed/tcsl248_processed.csv") 