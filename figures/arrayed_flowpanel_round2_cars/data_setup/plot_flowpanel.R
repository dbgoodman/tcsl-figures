#' Plot flow cytometry panel data
#' 
#' This script contains functions for plotting flow cytometry panel data
#' from TCSL248 and TCSL250 experiments.

library(data.table)
library(ggplot2)
library(patchwork)
library(here)

# Source the shared CAR color definitions
source(here("figures", "shared", "car_colors.R"))

#' Make a frequency plot for a specific subset and gate
#'
#' @param df The processed data.table
#' @param this_subset The subset to plot (e.g., "CD4", "CD8")
#' @param this_gate The gate to plot
#' @return A ggplot object
make_freq_plot <- function(df, this_subset, this_gate) {
  # Create the plot
  p <- ggplot(
    df[subset==this_subset & gate_clean %in% this_gate
      & xor(CAR == TRUE, car_name == "U")][, 
        car_stim := interaction(stimmed, car_name)][, 
        car_stim := reorder(car_stim, freq)]) + 
    geom_point(aes(x=car_stim, y=freq, color=rep_lbl)) + 
    # label TRUE as 'Stimmed' and FALSE as 'Unstimmed':
    facet_wrap(~stimmed, scales='free',
      labeller = as_labeller(c("TRUE" = 'Stimmed', "FALSE" = 'Unstimmed'))) +
    scale_x_discrete(labels = function(x) gsub("^[^\\.]*\\.", "", x)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = paste0(this_subset, ' - ', this_gate))
  
  return(p)
}

#' Make a count plot for a specific subset
#'
#' @param df The processed data.table
#' @param this_subset The subset to plot (e.g., "CD4", "CD8")
#' @return A ggplot object
make_count_plot <- function(df, this_subset) {
  # Create the plot
  p <- ggplot(
    df[subset=='CD3' & gate_clean == this_subset
      & xor(CAR == TRUE, car_name == "U")][, 
        car_rep := interaction(rep_type, car_name)][, 
        car_rep := reorder(car_rep, count_per_bead)]) + 
    geom_point(aes(x=car_rep, y=count_per_bead, color=rep_lbl)) + 
    facet_wrap(~rep_type, scales='free') +
    scale_x_discrete(labels = function(x) gsub("^[^\\.]*\\.", "", x)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = paste0(this_subset))
  
  return(p)
}

#' Make a mean MFI plot for a specific subset and gate
#'
#' @param df The processed data.table
#' @param this_subset The subset to plot (e.g., "CD4", "CD8")
#' @param this_gate The gate to plot
#' @param measure The measure to plot (default: "mean_mfi")
#' @return A ggplot object
make_mean_mfi_plot <- function(df, this_subset, this_gate, measure='mean_mfi') {
  # Create the plot
  p <- ggplot(
    df[subset==this_subset & gate_clean %in% this_gate
      & xor(CAR == TRUE, car_name == "U")][, 
        car_stim := interaction(stimmed, car_name)][, 
        car_stim := reorder(car_stim, get(measure))]) + 
    geom_point(aes(x=car_stim, y=get(measure), color=rep_lbl)) + 
    facet_wrap(~stimmed, scales='free') +
    scale_x_discrete(labels = function(x) gsub("^[^\\.]*\\.", "", x)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = paste0(this_subset, ' - ', this_gate))
  
  return(p)
}

#' Make a frequency plot for a specific subset and gate category
#'
#' @param df The processed data.table
#' @param this_subset The subset to plot (e.g., "CD4", "CD8")
#' @param this_category The gate category to plot (e.g., "Activation", "Differentiation", "Exhaustion")
#' @return A ggplot object
make_category_freq_plot <- function(df, this_subset, this_category, normalize = FALSE) {
  # Filter data for the specified subset and gate category
  plot_data <- df[subset == this_subset & gate_category == this_category]
  
  # Skip if no data
  if (nrow(plot_data) == 0) {
    message(sprintf("No data for %s %s, skipping", this_subset, this_category))
    return(NULL)
  }
  
  # Check if experiment column exists, if not, add a dummy one
  if (!"experiment" %in% names(plot_data)) {
    plot_data[, experiment := "experiment"]
  }
  
  # Normalize to zeta if requested
  if (normalize) {
    # For each experiment, subset, and gate_clean, calculate the zeta value
    zeta_values <- plot_data[car_name == "Zeta", 
                           .(zeta_value = mean(freq, na.rm = TRUE)), 
                           by = .(experiment, subset, gate_clean)]
    
    # Merge zeta values back to the data
    plot_data <- merge(plot_data, zeta_values, by = c("experiment", "subset", "gate_clean"), all.x = TRUE)
    
    # Calculate normalized values
    plot_data[, freq_norm := freq / zeta_value]
    
    # Use normalized values for plotting
    y_col <- "freq_norm"
    y_label <- "Frequency (normalized to Zeta)"
    title_suffix <- " (Normalized to Zeta)"
  } else {
    # Use raw frequency values
    y_col <- "freq"
    y_label <- "Frequency (%)"
    title_suffix <- ""
  }
  
  # Create the plot using the same approach as car_expression plots
  p <- ggplot(plot_data) +
    geom_point(aes(x = reorder(car_name, get(y_col), mean), 
                  y = get(y_col), 
                  color = experiment, 
                  shape = rep_lbl), 
              size = 3) +
    facet_wrap(~ gate_clean, scales = "free_y") +
    scale_y_continuous(name = y_label) +
    scale_color_manual(values = c("tcsl248" = "#E69F00", "tcsl250" = "#56B4E9"),
                      labels = c("tcsl248" = "TCSL248", "tcsl250" = "TCSL250")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "right") +
    labs(title = paste0(this_subset, " - ", this_category, title_suffix),
         x = "CAR Construct")
  
  if (normalize) {
    # Add horizontal line at y=1 for normalized plots
    p <- p + geom_hline(yintercept = 1, linetype = "dashed", color = "gray50")
  }
  
  return(p)
}

#' Make a mean MFI plot for a specific subset and gate category
#'
#' @param df The processed data.table
#' @param this_subset The subset to plot (e.g., "CD4", "CD8")
#' @param this_category The gate category to plot (e.g., "Exhaustion")
#' @param measure The measure to plot (default: "mean_mfi")
#' @return A ggplot object
make_category_mfi_plot <- function(df, this_subset, this_category, measure = "mean_mfi", normalize = FALSE) {
  # Filter data for the specified subset and gate category
  plot_data <- df[subset == this_subset & gate_category == this_category]
  
  # Skip if no data
  if (nrow(plot_data) == 0) {
    message(sprintf("No data for %s %s, skipping", this_subset, this_category))
    return(NULL)
  }
  
  # Check if experiment column exists, if not, add a dummy one
  if (!"experiment" %in% names(plot_data)) {
    plot_data[, experiment := "experiment"]
  }
  
  # Normalize to zeta if requested
  if (normalize) {
    # For each experiment, subset, and gate_clean, calculate the zeta value
    zeta_values <- plot_data[car_name == "Zeta", 
                           .(zeta_value = mean(get(measure), na.rm = TRUE)), 
                           by = .(experiment, subset, gate_clean)]
    
    # Merge zeta values back to the data
    plot_data <- merge(plot_data, zeta_values, by = c("experiment", "subset", "gate_clean"), all.x = TRUE)
    
    # Calculate normalized values
    plot_data[, mfi_norm := get(measure) / zeta_value]
    
    # Use normalized values for plotting
    y_col <- "mfi_norm"
    y_label <- "MFI (normalized to Zeta)"
    title_suffix <- " (Normalized to Zeta)"
  } else {
    # Use raw MFI values
    y_col <- measure
    y_label <- "Mean Fluorescence Intensity"
    title_suffix <- ""
  }
  
  # Create the plot using the same approach as car_expression plots
  p <- ggplot(plot_data) +
    geom_point(aes(x = reorder(car_name, get(y_col), mean), 
                  y = get(y_col), 
                  color = experiment, 
                  shape = rep_lbl), 
              size = 3) +
    facet_wrap(~ gate_clean, scales = "free_y") +
    scale_y_continuous(name = y_label) +
    scale_color_manual(values = c("tcsl248" = "#E69F00", "tcsl250" = "#56B4E9"),
                      labels = c("tcsl248" = "TCSL248", "tcsl250" = "TCSL250")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "right") +
    labs(title = paste0(this_subset, " - ", this_category, " MFI", title_suffix),
         x = "CAR Construct")
  
  if (normalize) {
    # Add horizontal line at y=1 for normalized plots
    p <- p + geom_hline(yintercept = 1, linetype = "dashed", color = "gray50")
  }
  
  return(p)
}

#' Plot chronic stimulation expansion after 7 days
#'
#' @param df The processed data.table
#' @return A ggplot object
plot_chronic_stim_expansion <- function(df) {
  make_count_plot(df[rep != 'X' & car_name != 'U'], 'CD8') /
  make_count_plot(df[rep != 'X' & car_name != 'U'], 'CD4')
}

#' Plot percent EGFR+
#'
#' @param df The processed data.table
#' @return A ggplot object
plot_percent_egfr <- function(df) {
  make_freq_plot(df, 'CD3', 'live_car')
}

#' Plot bead-based counts
#'
#' @param df The processed data.table
#' @param bead_count_file Path to the bead count file
#' @return A list of ggplot objects
plot_bead_based_counts <- function(df, bead_count_file) {
  # Extract flow counts
  flow_counts <- df[rep != 'X' & car_name != 'U', 
    .(rs_bead_count = mean(
      count_per_bead[gate_clean=='live_untr'] + count_per_bead[gate_clean=='live_car'])),
    by=.(car_id, well_id, rep, rep_lbl, car_name)]
  
  # Load bead count data
  bead_count_df <- fread(bead_count_file)
  if ("V1" %in% names(bead_count_df)) {
    bead_count_df[, V1 := NULL]
  }
  setnames(bead_count_df, c('name'), c('well_id'))
  
  # Merge with flow counts
  bead_count_df <- bead_count_df[flow_counts, on='well_id']
  
  # Create plots
  plots <- list()
  
  # Plot 1: Raw values by variable
  plots$raw_by_variable <- ggplot(
    melt(bead_count_df, id.vars=c('well_id', 'car_id','rep', 'car_name','rep_lbl'))[
      !grepl('std',variable)][, car_rep := reorder(interaction(car_name, rep_lbl, variable, sep='|'), value)]) + 
    geom_point(aes(x=car_rep, y=value, color=rep_lbl)) + 
    facet_wrap(rep_lbl~variable, scales='free') + 
    scale_x_discrete(labels = function(x) gsub("^([^\\|]+)\\|.*", "\\1", x)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  # Create melted data frame with scaled values
  bead_count_df_melt <- melt(bead_count_df, id.vars=c('well_id', 'car_id','rep', 'rep_lbl', 'car_name'))[
    !grepl('std',variable)][, 
      scaled_value := scale(value), by=.(variable)][, 
      car_rep := reorder(interaction(car_name, rep_lbl, sep='|'), scaled_value)]
  
  # Plot 2: Scaled values by replicate
  plots$scaled_by_replicate <- ggplot(bead_count_df_melt) + 
    geom_point(aes(x=car_rep, y=scaled_value, color=rep_lbl)) + 
    facet_wrap(~rep_lbl, scales='free') + 
    scale_x_discrete(labels = function(x) gsub("^([^\\|]+)\\|.*", "\\1", x)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  # Plot 3: Scaled values by CAR
  plots$scaled_by_car <- ggplot(bead_count_df_melt[, car_name := reorder(car_name, scaled_value)][, scaled_value_mean := mean(scaled_value), by=.(car_name)]) + 
    geom_point(aes(x=car_name, y=scaled_value, color=rep_lbl)) + 
    geom_point(aes(x=car_name, y=scaled_value_mean), size=3) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x='CAR', y='Scaled Total Count')
  
  # Calculate CAR percentages
  car_pcts <- df[
    rep != 'X' & car_name != 'U' & subset == 'CD3' & gate_clean %in% c('live_car','live_untr'), 
    .(pct_car= mean(
      freq[gate_clean == 'live_car']/(freq[gate_clean == 'live_car'] + 
        freq[gate_clean == 'live_untr']))), 
    by=.(car_id, well_id, rep, car_name)]
  
  # Add CAR percentages to melted data
  bead_count_df_melt <- bead_count_df_melt[
    car_pcts, on=.(car_id, well_id, rep, car_name)][, value_car_pos := value * pct_car][, 
        scaled_value_car_pos := scale(value_car_pos), by=.(variable)]
  
  # Plot 4: Scaled CAR+ values by replicate
  plots$scaled_car_pos_by_replicate <- ggplot(bead_count_df_melt[, car_rep := reorder(car_rep, scaled_value_car_pos)]) + 
    geom_point(aes(x=car_rep, y=scaled_value_car_pos, color=rep_lbl)) + 
    facet_wrap(~rep_lbl, scales='free') + 
    scale_x_discrete(labels = function(x) gsub("^([^\\|]+)\\|.*", "\\1", x)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x='CAR', y='Scaled CAR+ Count')
  
  # Plot 5: Scaled CAR+ values by CAR
  plots$scaled_car_pos_by_car <- ggplot(bead_count_df_melt[, car_name := reorder(car_name, scaled_value_car_pos)][, 
      scaled_value_car_pos_mean := mean(scaled_value_car_pos), by=.(car_name)]) + 
    geom_point(aes(x=car_name, y=scaled_value_car_pos, color=rep_lbl)) + 
    geom_point(aes(x=car_name, y=scaled_value_car_pos_mean), size=3) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x='CAR', y='Scaled CAR+ Count')
  
  return(plots)
}

#' Plot CD8 differentiation by category
#'
#' @param df The processed data.table
#' @return A ggplot object
plot_cd8_differentiation_by_category <- function(df) {
  make_category_freq_plot(df, 'CD8', 'Differentiation')
}

#' Plot CD8 exhaustion frequency by category
#'
#' @param df The processed data.table
#' @return A ggplot object
plot_cd8_exhaustion_freq_by_category <- function(df) {
  make_category_freq_plot(df, 'CD8', 'Exhaustion')
}

#' Plot CD8 exhaustion MFI by category
#'
#' @param df The processed data.table
#' @return A ggplot object
plot_cd8_exhaustion_mfi_by_category <- function(df) {
  make_category_mfi_plot(df, 'CD8', 'Exhaustion')
}

#' Plot CD8 activation by category
#'
#' @param df The processed data.table
#' @return A ggplot object
plot_cd8_activation_by_category <- function(df) {
  make_category_freq_plot(df, 'CD8', 'Activation')
}

#' Plot CD4 differentiation by category
#'
#' @param df The processed data.table
#' @return A ggplot object
plot_cd4_differentiation_by_category <- function(df) {
  make_category_freq_plot(df, 'CD4', 'Differentiation')
}

#' Plot CD4 exhaustion frequency by category
#'
#' @param df The processed data.table
#' @return A ggplot object
plot_cd4_exhaustion_freq_by_category <- function(df) {
  make_category_freq_plot(df, 'CD4', 'Exhaustion')
}

#' Plot CD4 exhaustion MFI by category
#'
#' @param df The processed data.table
#' @return A ggplot object
plot_cd4_exhaustion_mfi_by_category <- function(df) {
  make_category_mfi_plot(df, 'CD4', 'Exhaustion')
}

#' Plot CD4 activation by category
#'
#' @param df The processed data.table
#' @return A ggplot object
plot_cd4_activation_by_category <- function(df) {
  make_category_freq_plot(df, 'CD4', 'Activation')
}

#' Plot CD3 differentiation by category
#'
#' @param df The processed data.table
#' @return A ggplot object
plot_cd3_differentiation_by_category <- function(df) {
  make_category_freq_plot(df, 'CD3', 'Differentiation')
}

#' Plot CD3 exhaustion frequency by category
#'
#' @param df The processed data.table
#' @return A ggplot object
plot_cd3_exhaustion_freq_by_category <- function(df) {
  make_category_freq_plot(df, 'CD3', 'Exhaustion')
}

#' Plot CD3 exhaustion MFI by category
#'
#' @param df The processed data.table
#' @return A ggplot object
plot_cd3_exhaustion_mfi_by_category <- function(df) {
  make_category_mfi_plot(df, 'CD3', 'Exhaustion')
}

#' Plot CD3 activation by category
#'
#' @param df The processed data.table
#' @return A ggplot object
plot_cd3_activation_by_category <- function(df) {
  make_category_freq_plot(df, 'CD3', 'Activation')
}

#' Create a tile plot of marker expression
#'
#' @param df The processed data.table
#' @param subset_filter Optional filter for subset (e.g., "CD4", "CD8", "CD3")
#' @return A list of ggplot objects
create_tile_plots <- function(df, subset_filter = NULL) {
  # Define gates and categories
  gates <- c('CD25','CD27','CD127','Eff','Emra','TIM3','Mem','Naive','CD39','PD1','LAG3')
  
  # Ensure gate categories are assigned
  df[gate_clean %in% c('CD25','CD27','CD127'), gate_category := 'Activation']
  df[gate_clean %in% c('Eff','Emra','Mem','Naive'), gate_category := 'Differentiation']
  df[gate_clean %in% c('TIM3','LAG3','CD39','PD1'), gate_category := 'Exhaustion']
  
  # Filter by subset if specified
  if (!is.null(subset_filter)) {
    df_subset <- df[subset == subset_filter]
  } else {
    df_subset <- df
  }
  
  # Create plots
  plots <- list()
  
  # Frequency tile plot
  freq_data <- df_subset[car_name != 'U' & subset %in% c('CD4','CD8','CD3') & gate_clean %in% gates
    & xor(CAR == TRUE, car_name == "U"), .(mean_freq= mean(freq)), 
    by=.(car_name, subset, gate_clean, gate_category)][,
          scaled_mean_freq := scale(mean_freq), by=.(subset, gate_clean)][, 
            car_name_cat := reorder(interaction(car_name, gate_category, subset, sep='|'), scaled_mean_freq)]
  
  # Create the frequency tile plot
  p_freq <- ggplot(freq_data) + 
    geom_tile(aes(x=car_name_cat, fill=scaled_mean_freq, y=gate_clean)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_x_discrete(labels = function(x) gsub("^([^\\|]*)\\|.*", "\\1", x)) + 
    scale_fill_distiller(palette='Spectral') +
    geom_text(aes(x=car_name_cat, y=gate_clean, label=round(mean_freq))) + 
    facet_wrap(subset~gate_category, scales='free')
  
  plots$freq_tile <- p_freq
  
  # MFI tile plot for exhaustion markers
  mfi_data <- df_subset[car_name != 'U' & subset %in% c('CD4','CD8','CD3') & gate_category %in% 'Exhaustion'
    & xor(CAR == TRUE, car_name == "U"), .(mean_mfi= mean(mean_mfi)), 
    by=.(car_name, subset, gate_clean, gate_category)][,
          scaled_mean_mfi := scale(mean_mfi), by=.(subset, gate_clean)][, 
            car_name_cat := reorder(interaction(car_name, gate_category, subset, sep='|'), scaled_mean_mfi)]
  
  # Create the MFI tile plot
  p_mfi <- ggplot(mfi_data) + 
    geom_tile(aes(x=car_name_cat, fill=scaled_mean_mfi, y=gate_clean)) +
    scale_x_discrete(labels = function(x) gsub("^([^\\|]*)\\|.*", "\\1", x)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_distiller(palette='Spectral') +
    geom_text(aes(x=car_name_cat, y=gate_clean, label=round(mean_mfi))) + 
    facet_wrap(subset~gate_category, scales='free')
  
  plots$mfi_tile <- p_mfi
  
  return(plots)
}

#' Save all plots for a specific experiment
#'
#' @param df The processed data.table
#' @param experiment_id The experiment ID (e.g., "tcsl248" or "tcsl250")
#' @param bead_count_file Path to the bead count file (optional)
#' @param output_dir Directory to save plots (default: figures/arrayed_flowpanel_round2_cars/plots/{experiment_id})
#' @return NULL
save_all_plots <- function(df, experiment_id, bead_count_file = NULL, output_dir = NULL) {
  # Set output directory if not provided
  if (is.null(output_dir)) {
    output_dir <- here("figures", "arrayed_flowpanel_round2_cars", "plots", experiment_id)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  print(unique(df[, .(car_name, CAR, stimmed, subset)]))
  
  # Ensure gate categories are assigned
  df[gate_clean %in% c('CD25','CD27','CD127'), gate_category := 'Activation']
  df[gate_clean %in% c('Eff','Emra','Mem','Naive'), gate_category := 'Differentiation']
  df[gate_clean %in% c('TIM3','LAG3','CD39','PD1'), gate_category := 'Exhaustion']
  
  # Save chronic stim expansion plot
  ggsave(file.path(output_dir, "chronic_stim_expansion.pdf"), 
         plot_chronic_stim_expansion(df), 
         width = 10, height = 10)
  
  # Save percent EGFR+ plot
  ggsave(file.path(output_dir, "percent_egfr.pdf"), 
         plot_percent_egfr(df), 
         width = 10, height = 6)
  
  # Save bead-based count plots if bead count file is provided
  if (!is.null(bead_count_file)) {
    bead_plots <- plot_bead_based_counts(df, bead_count_file)
    for (name in names(bead_plots)) {
      ggsave(file.path(output_dir, paste0("bead_counts_", name, ".pdf")), 
             bead_plots[[name]], 
             width = 10, height = 8)
    }
  }
  
  # Save CD8 plots by category
  ggsave(file.path(output_dir, "cd8_differentiation_by_category.pdf"), 
         plot_cd8_differentiation_by_category(df), 
         width = 12, height = 10)
  
  ggsave(file.path(output_dir, "cd8_exhaustion_freq_by_category.pdf"), 
         plot_cd8_exhaustion_freq_by_category(df), 
         width = 12, height = 10)
  
  ggsave(file.path(output_dir, "cd8_exhaustion_mfi_by_category.pdf"), 
         plot_cd8_exhaustion_mfi_by_category(df), 
         width = 12, height = 10)
  
  ggsave(file.path(output_dir, "cd8_activation_by_category.pdf"), 
         plot_cd8_activation_by_category(df), 
         width = 12, height = 10)
  
  # Save CD4 plots by category
  ggsave(file.path(output_dir, "cd4_differentiation_by_category.pdf"), 
         plot_cd4_differentiation_by_category(df), 
         width = 12, height = 10)
  
  ggsave(file.path(output_dir, "cd4_exhaustion_freq_by_category.pdf"), 
         plot_cd4_exhaustion_freq_by_category(df), 
         width = 12, height = 10)
  
  ggsave(file.path(output_dir, "cd4_exhaustion_mfi_by_category.pdf"), 
         plot_cd4_exhaustion_mfi_by_category(df), 
         width = 12, height = 10)
  
  ggsave(file.path(output_dir, "cd4_activation_by_category.pdf"), 
         plot_cd4_activation_by_category(df), 
         width = 12, height = 10)
  
  # Save CD3 plots by category
  ggsave(file.path(output_dir, "cd3_differentiation_by_category.pdf"), 
         plot_cd3_differentiation_by_category(df), 
         width = 12, height = 10)
  
  ggsave(file.path(output_dir, "cd3_exhaustion_freq_by_category.pdf"), 
         plot_cd3_exhaustion_freq_by_category(df), 
         width = 12, height = 10)
  
  ggsave(file.path(output_dir, "cd3_exhaustion_mfi_by_category.pdf"), 
         plot_cd3_exhaustion_mfi_by_category(df), 
         width = 12, height = 10)
  
  ggsave(file.path(output_dir, "cd3_activation_by_category.pdf"), 
         plot_cd3_activation_by_category(df), 
         width = 12, height = 10)
  
  # Save tile plots for all subsets
  tile_plots <- create_tile_plots(df)
  for (name in names(tile_plots)) {
    ggsave(file.path(output_dir, paste0("tile_plot_", name, ".pdf")), 
           tile_plots[[name]], 
           width = 20, height = 8)
  }
  
  # Save CD3-only tile plots
  cd3_tile_plots <- create_tile_plots(df, subset_filter = "CD3")
  for (name in names(cd3_tile_plots)) {
    ggsave(file.path(output_dir, paste0("cd3_tile_plot_", name, ".pdf")), 
           cd3_tile_plots[[name]], 
           width = 20, height = 8)
  }
  
  message(sprintf("All plots saved to %s", output_dir))
} 