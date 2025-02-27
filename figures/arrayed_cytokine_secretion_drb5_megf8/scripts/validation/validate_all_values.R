library(data.table)
library(here)

# Function to process data exactly as in the Rmd
process_rmd_style <- function(data_dir, populations_file) {
  df_raw <- fread(file.path(data_dir, populations_file))
  
  # keep *.fcs only
  df_raw <- df_raw[grepl('.fcs', V1)]
  df_raw[, `:=`(plate = as.integer(sub(".*Plate_([0-9]+).*", "\\1", V1)),
            well = sub(".*-([A-Z][0-9]+) .*", "\\1", V1))]
  df_raw[, grep("^V", names(df_raw)) := NULL]
  df_melt <- unique(melt(df_raw, id.vars=c('plate', 'well'), variable.name='gate_str', value.name='value'))
  df_melt[, `:=`(
    gate = sub(" \\|.*", "", sub(".*/", "", gate_str)),  
    parent_gate = fifelse(grepl("^([^/]+) \\|", gate_str),        
                   NA_character_,                         
                   gsub(".*/([^/]+)/[^/]+ \\|.*", "\\1", gsub('Q10: IFNg\\+ , TNF-a\\+/', '', gate_str))),  
    measure = ifelse(grepl("\\| Freq", gate_str), "freq", 
                    ifelse(grepl("\\| Count", gate_str), "count", NA))
  )]
  
  df_melt[parent_gate == 'Single Cells', parent_gate := 'Lymphocytes']
  df_melt <- df_melt[is.na(parent_gate) | parent_gate != 'Lymphocytes' | grepl('Freq. of Lymphocytes', gate_str)]
  df_melt[, CAR := grepl('_car', gate_str)]
  df_melt[, subset := 'None']
  df_melt[grepl('/(CD[348])/', gate_str), subset := gsub('.*/(CD[348])/.*', '\\1', gate_str)]
  df_melt[, parent_subset := subset]
  df_melt[grepl('^CD[348]$', parent_gate), parent_subset := 'None']
  df_melt[, parent_CAR := CAR]
  df_melt[gate == 'live_car', parent_CAR := FALSE]
  df_melt <- dcast(df_melt, ... ~ measure, value.var='value')
  df_melt[, gate_str := sub(' \\|.*', '', gate_str)]
  df_melt[is.na(parent_gate), freq := freq[!is.na(freq)], by= .(plate, well, gate, CAR, subset)]
  df_melt[is.na(parent_gate), count := count[!is.na(count)], by= .(plate, well, gate, CAR, subset)]
  df_melt <- unique(df_melt)
  df_melt[is.na(freq), freq := 0]
  
  while (any(is.na(df_melt[!(gate %in% c('beads','Lymphocytes')), count]))) {
    df_new_counts <- df_melt[!is.na(count), .(CAR, gate, subset, well, plate, count)][df_melt[is.na(count)],   
      on=.(CAR=parent_CAR, gate=parent_gate, subset=parent_subset, well, plate)][
        !is.na(count)][, 
        .(plate, well, gate=i.gate, CAR=i.CAR, subset=i.subset,
          count=round(freq/100 * count))]
    
    df_melt[df_new_counts, 
            on = .(plate, well, gate, CAR, subset), 
            count := fifelse(is.na(count), i.count, count)]
  }
  
  df_melt[, gate_clean := gate]
  df_melt[, gate_clean := gsub("TNF-a", "TNFa", gate_clean)]
  df_melt[, gate_clean := gsub("\\+", "_pos", gate_clean)]
  df_melt[, gate_clean := gsub("\\-", "_neg", gate_clean)]
  df_melt[, gate_clean := gsub(" , ", ".", gate_clean)]
  df_melt[, gate_clean := gsub("^\\s+", "", gate_clean)]
  df_melt[, gate_clean := gsub("^Q[0-9]+: ", "", gate_clean)]
  df_melt[grepl('Q10: IFNg\\+ , TNF-a\\+/', gate_str), 
    gate_clean := paste0('IFNg_pos.TNFa_pos.', gate_clean)]
  
  # Load plate map
  plate_map <- fread(file.path(data_dir, "plate_map.csv"))
  plate_map[, rn := LETTERS[1:nrow(plate_map)]]
  long_plate <- melt(plate_map, id.vars = "rn", variable.name = "col", value.name = "car_id")
  long_plate[, `:=`(well = paste0(rn, as.integer(sub("V", "", col))),
                    row = rn,
                    col = as.integer(sub("V", "", col)))]
  long_plate <- long_plate[!is.na(car_id)]
  df_melt <- merge(df_melt, long_plate[, .(well, car_id, row, col)], by = "well", all.x = TRUE)
  df_melt[, car_id := factor(car_id, levels = c(as.character(1:16), "NEG","POS"))]
  
  # Load CAR map
  car_map <- fread(file.path(data_dir, "car_map.csv"))
  car_map[, car_id := factor(car_id)]
  car_map[, car_name := factor(car_name)]
  df_melt <- merge(df_melt, car_map, by = "car_id", all.x = TRUE)
  df_melt[is.na(car_name), car_name := car_id]
  
  return(df_melt)
}

# Process data using Rmd style
tcsl248_rmd <- process_rmd_style('/Users/dbg/Library/CloudStorage/Box-Box/tcsl/flow/20241112_tcsl248_cytokine', 'populations.csv')
tcsl250_rmd <- process_rmd_style('/Users/dbg/Library/CloudStorage/Box-Box/tcsl/flow/20241213.tcsl250 cytokines', 'populations.csv')

# Load our processed data
tcsl248_new <- fread(here('data', 'processed', 'cytokine_secretion', 'tcsl248_processed.csv'))
tcsl250_new <- fread(here('data', 'processed', 'cytokine_secretion', 'tcsl250_processed.csv'))

# Function to compare datasets
compare_datasets <- function(df1, df2, name1, name2) {
  cat("\nComparing", name1, "vs", name2, "\n")
  cat("Dimensions:", dim(df1), "vs", dim(df2), "\n")
  
  # Compare column names
  col_diff1 <- setdiff(names(df1), names(df2))
  col_diff2 <- setdiff(names(df2), names(df1))
  if (length(col_diff1) > 0) cat("Columns in", name1, "but not in", name2, ":", col_diff1, "\n")
  if (length(col_diff2) > 0) cat("Columns in", name2, "but not in", name1, ":", col_diff2, "\n")
  
  # Compare values for key columns
  key_cols <- c("car_name", "gate_clean", "subset", "freq", "count")
  for (col in intersect(key_cols, names(df1))) {
    # Sort both datasets to ensure matching order
    df1_sorted <- df1[order(car_name, gate_clean, subset, well)]
    df2_sorted <- df2[order(car_name, gate_clean, subset, well)]
    
    # Compare values
    if (is.numeric(df1[[col]])) {
      max_diff <- max(abs(df1_sorted[[col]] - df2_sorted[[col]]), na.rm=TRUE)
      cat("\nMax difference in", col, ":", max_diff, "\n")
      if (max_diff > 0.01) {
        # Show examples of differences
        diffs <- which(abs(df1_sorted[[col]] - df2_sorted[[col]]) > 0.01)
        if (length(diffs) > 0) {
          cat("Example differences in", col, ":\n")
          print(data.table(
            car_name = df1_sorted$car_name[diffs[1:min(5, length(diffs))]],
            gate_clean = df1_sorted$gate_clean[diffs[1:min(5, length(diffs))]],
            subset = df1_sorted$subset[diffs[1:min(5, length(diffs))]],
            val1 = df1_sorted[[col]][diffs[1:min(5, length(diffs))]],
            val2 = df2_sorted[[col]][diffs[1:min(5, length(diffs))]]
          ))
        }
      }
    } else {
      # For non-numeric columns, compare unique values
      uniq1 <- sort(unique(df1_sorted[[col]]))
      uniq2 <- sort(unique(df2_sorted[[col]]))
      if (!identical(uniq1, uniq2)) {
        cat("\nDifferences in unique values for", col, ":\n")
        cat("In", name1, "but not in", name2, ":", setdiff(uniq1, uniq2), "\n")
        cat("In", name2, "but not in", name1, ":", setdiff(uniq2, uniq1), "\n")
      }
    }
  }
}

# Compare TCSL248 datasets
compare_datasets(tcsl248_rmd, tcsl248_new, "TCSL248 RMD", "TCSL248 New")

# Compare TCSL250 datasets
compare_datasets(tcsl250_rmd, tcsl250_new, "TCSL250 RMD", "TCSL250 New")

# Compare summary statistics for plotting
chosen_gate <- c(
  "IL2_pos", "TNFa_pos", "IFNg_pos",
  "IFNg_pos.IL2_pos", "IFNg_pos.TNFa_pos",
  "IL2_pos.TNFa_pos", "IFNg_pos.TNFa_pos.IL2_pos")

get_plot_summary <- function(df) {
  df[CAR == TRUE & !is.na(car_name)][
    gate_clean %in% chosen_gate, .(mean_freq = mean(freq)), 
    by=.(car_name, gate_clean)]
}

cat("\nComparing plot summary statistics:\n")
summary_248_rmd <- get_plot_summary(tcsl248_rmd)
summary_248_new <- get_plot_summary(tcsl248_new)
summary_250_rmd <- get_plot_summary(tcsl250_rmd)
summary_250_new <- get_plot_summary(tcsl250_new)

cat("\nTCSL248 plot values comparison:\n")
print(merge(summary_248_rmd, summary_248_new, 
           by=c("car_name", "gate_clean"), 
           suffixes=c("_rmd", "_new"))[abs(mean_freq_rmd - mean_freq_new) > 0.01])

cat("\nTCSL250 plot values comparison:\n")
print(merge(summary_250_rmd, summary_250_new, 
           by=c("car_name", "gate_clean"), 
           suffixes=c("_rmd", "_new"))[abs(mean_freq_rmd - mean_freq_new) > 0.01]) 