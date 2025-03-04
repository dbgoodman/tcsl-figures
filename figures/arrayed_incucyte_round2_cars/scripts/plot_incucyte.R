#' Plot incucyte data for TCSL248 and TCSL250 experiments
#'
#' This script creates visualizations of the incucyte data for the arrayed CAR experiments.

library(data.table)
library(ggplot2)
library(patchwork)
library(stringr)
library(ggrepel)
library(viridis)
library(here)
library(fs)

# Source the shared CAR color definitions
source(here("figures", "shared", "car_colors.R"))

# Define paths
processed_dir <- here("data", "processed", "incucyte_round2")
figure_dir <- here("figures", "arrayed_incucyte_round2_cars")

# Load processed data
all_data <- fread(file.path(processed_dir, "all_incucyte_data.csv"))

# Set theme for all plots
theme_set(theme_bw() +
            theme(
              panel.grid.minor = element_blank(),
              strip.background = element_rect(fill = "white"),
              strip.text = element_text(face = "bold"),
              legend.position = "bottom"
            ))

# Define custom color palette for CAR variants
# Group similar CARs with related colors
# car_colors <- c(
#   # Zeta - Black
#   "Zeta" = "#000000",
#   
#   # CD28 variants - Reds
#   "CD28" = "#E41A1C",
#   "CD28.lev11" = "#FB6A4A",
#   "CD28.lev20" = "#FCAE91",
#   
#   # TNR9 variants - Greens
#   "TNR9" = "#33A02C",
#   "TNR9.lev8" = "#78C679",
#   "TNR9.lev16" = "#C2E699",
#   
#   # DRB5/DRB528C variants - Blues
#   "DRB5" = "#0868AC",
#   "DRB5.lev8" = "#43A2CA",
#   "DRB528C" = "#1F78B4",
#   "DRB528C.lev26" = "#7BCCC4",
#   
#   # MEGF8 variants - Oranges
#   "MEGF8" = "#FF7F00",
#   "MEGF8.lev30" = "#FDB462",
#   
#   # Other CARs - Various colors
#   "MBTP1" = "#FFFF33",
#   "COMB.opt.16aa" = "#A65628",
#   "COMB.opt.27aa" = "#F781BF",
#   
#   # Controls - Gray
#   "U" = "#999999",
#   "K" = "#CCCCCC"
# )

# Define final timepoints for each experiment
final_timepoints_map <- list(
  tcsl248 = list(
    "week_0" = 60,
    "week_1" = 80
  ),
  tcsl250 = list(
    "week_0" = 60,
    "week_1" = 80,
    "week_0_purified" = 60
  )
)


# Modify create_time_course_plot function to handle combined experiments
create_time_course_plot <- function(data, title, facet_ncol = 4, car_subset = NULL) {
  # Filter data to just what we need
  plot_data <- data[measurement == "incucyte_total_int_intensity"]
  
  # Apply CAR subset if provided
  if (!is.null(car_subset)) {
    plot_data <- plot_data[car_name %in% car_subset]
  }
  
  # Filter out control wells
  plot_data <- plot_data[!(car_name %in% c('K','U'))]
  
  # Create a combined dataset with appropriate filtering by final timepoints
  filtered_data <- data.table()

  # Process each experiment-week combination and filter by appropriate final timepoint
  for (exp in unique(plot_data$experiment)) {
    exp_lower <- tolower(exp)
    for (wk in unique(plot_data[experiment == exp]$week)) {
      # Get final timepoint from the map
      week_key <- ifelse(wk == "0_purified", "week_0_purified", paste0("week_", wk))
      final_timepoint <- final_timepoints_map[[exp_lower]][[week_key]]
      
      # Filter data for this timepoint
      timepoint_data <- plot_data[experiment == exp & week == wk & Elapsed <= final_timepoint]
      
      # Create a combined facet label with experiment, week, ratio, and final timepoint
      timepoint_data[, facet_label := paste(
        exp, 
        ifelse(wk == "0_purified", "Week 0 Purified", paste("Week", wk)),
        paste0(ratio, " A549 cells : 1 CAR-T cell"),
        paste0("Final: ", final_timepoint, " hours"),
        sep = "\n"
      )]
      
      # Add final timepoint as a column
      timepoint_data[, final_timepoint := final_timepoint]
      
      # Add week as numeric for sorting
      timepoint_data[, week_num := as.numeric(gsub("_purified", "", wk))]
      
      # Add to the filtered dataset
      filtered_data <- rbind(filtered_data, timepoint_data)
    }
  }
  
  # Order the facet labels based on week, experiment, and ratio
  filtered_data[, facet_label := factor(
    facet_label,
    levels = unique(filtered_data[order(week_num, experiment, as.numeric(ratio))]$facet_label)
  )]
  
  # Create the plot
  p <- ggplot(filtered_data) + 
    geom_line(aes(x = Elapsed, y = mean_value_none, 
                  group = car_name, 
                  color = car_name)) + 
    facet_wrap(~ facet_label, scales = "free", ncol = facet_ncol) +
    theme_minimal() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      legend.position = "bottom",
      panel.spacing = unit(1, "lines"),
      strip.text = element_text(face = "bold")
    ) +
    labs(title = title,
         x = 'Elapsed Time (hours)', 
         y = 'Fraction Cells Remaining\n(via normalized intensity vs no T cells)',
         color = "CAR") +
    scale_color_manual(values = car_colors) +
    # Extend x-axis to allow space for labels (20% extra on right side)
    scale_x_continuous(
      expand = expansion(mult = c(0.02, 0.2))  # Expand left by 2%, right by 20%
    )
  
  # Add vertical lines for final timepoints
  # We'll do this by creating a separate data frame for each facet
  vline_data <- unique(filtered_data[, .(facet_label, final_timepoint)])
  
  p <- p + geom_vline(data = vline_data, 
                     aes(xintercept = final_timepoint),
                     linetype = "dashed", color = "black", alpha = 0.5)
  
  # Get the last points for each car_name and facet
  last_points <- filtered_data[, .SD[which.max(Elapsed)], by = .(car_name, facet_label)]
  
  # Add labels
  p <- p + geom_text_repel(
    data = last_points,
    aes(x = Elapsed, y = mean_value_none, 
        label = car_name, 
        color = car_name,
        segment.color = after_scale(color)),
    hjust = 0,
    nudge_x = 5,
    direction = "y",
    segment.size = 0.5,
    segment.linetype = "dotted",
    box.padding = 0.2,
    point.padding = 0.1,
    force = 1,
    max.overlaps = 20,
    size = 3,
    min.segment.length = 0,
    seed = 42
  )
  
  return(p)
}

# Add a new function to create dot plots of final timepoint values
create_final_timepoint_dotplot <- function(data, title, car_subset = NULL) {
  # Filter data to just what we need
  plot_data <- data[measurement == "incucyte_total_int_intensity"]
  
  # Apply CAR subset if provided
  if (!is.null(car_subset)) {
    plot_data <- plot_data[car_name %in% car_subset]
  }
  
  # Filter out control wells
  plot_data <- plot_data[!(car_name %in% c('K','U'))]
  
  # Create a dataset with only the final timepoint values
  final_data <- data.table()
  
  # Process each experiment-week combination and get the final timepoint value
  for (exp in unique(plot_data$experiment)) {
    # Handle the renamed experiment for final timepoint lookup
    exp_lower <- tolower(gsub("_purified$", "", exp))
    
    for (wk in unique(plot_data[experiment == exp]$week)) {
      # Get final timepoint from the map
      # For tcsl250_purified, we need to use the purified timepoint
      if (exp == "tcsl250_purified") {
        week_key <- "week_0_purified"
      } else {
        week_key <- paste0("week_", wk)
      }
      
      final_timepoint <- final_timepoints_map[[exp_lower]][[week_key]]
      
      # Get data closest to the final timepoint for each CAR and ratio
      for (r in unique(plot_data[experiment == exp & week == wk]$ratio)) {
        for (car in unique(plot_data[experiment == exp & week == wk & ratio == r]$car_name)) {
          # Get data for this specific combination
          car_data <- plot_data[experiment == exp & week == wk & ratio == r & car_name == car]
          
          # Find the closest timepoint to the final timepoint
          closest_idx <- which.min(abs(car_data$Elapsed - final_timepoint))
          if (length(closest_idx) > 0) {
            final_point <- car_data[closest_idx]
            # Add final timepoint to the data
            final_point[, final_timepoint := final_timepoint]
            final_data <- rbind(final_data, final_point)
          }
        }
      }
    }
  }
  
  # Create a facet label for ratio and week with final timepoint
  final_data[, facet_label := paste(
    paste0(ratio, " A549 cells : 1 CAR-T cell"),
    paste("Week", week),
    paste0("Final: ", final_timepoint, " hours"),
    sep = "\n"
  )]
  
  # Add week as numeric for sorting
  final_data[, week_num := as.numeric(week)]
  
  # Order the facet labels by week then ratio
  final_data[, facet_label := factor(
    facet_label,
    levels = unique(final_data[order(week_num, as.numeric(ratio))]$facet_label)
  )]
  
  # Create the dot plot with experiment on x-axis
  p <- ggplot(final_data, aes(x = experiment, y = mean_value_none, color = car_name)) +
    geom_point(size = 3) +
    facet_wrap(~ facet_label, scales = "free", nrow = 1) +
    scale_color_manual(values = car_colors) +
    labs(
      title = title,
      x = "Experiment",
      y = "Fraction Cells Remaining at Final Timepoint",
      color = "CAR"
    ) +
    theme_minimal() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      legend.position = "bottom",
      panel.spacing = unit(1, "lines"),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  return(p)
}

# Initialize plots list
plots <- list()

# Define experiment-specific plots with their dimensions
experiment_plots <- list(
  tcsl248_wk0 = list(
    title = "TCSL248 Week 0 Incucyte Killing",
    dimensions = c(8, 12),
    data_filter = expression(experiment == "tcsl248" & week == 0),
    facet_ncol = 1
  ),
  tcsl248_wk1 = list(
    title = "TCSL248 Week 1 Incucyte Killing",
    dimensions = c(8, 12),
    data_filter = expression(experiment == "tcsl248" & week == 1),
    facet_ncol = 1
  ),
  tcsl250_wk0_standard = list(
    title = "TCSL250 Week 0 Standard Incucyte Killing",
    dimensions = c(8, 12),
    data_filter = expression(experiment == "tcsl250" & week == 0),
    facet_ncol = 1
  ),
  tcsl250_wk0_purified = list(
    title = "TCSL250 Week 0 Purified Incucyte Killing",
    dimensions = c(8, 8),
    data_filter = expression(experiment == "tcsl250" & week == "0_purified"),
    facet_ncol = 1
  ),
  tcsl250_wk1 = list(
    title = "TCSL250 Week 1 Incucyte Killing",
    dimensions = c(8, 6),
    data_filter = expression(experiment == "tcsl250" & week == 1),
    facet_ncol = 1
  )
)

# Define the CAR family groups

cars_to_compare <- c("Zeta", "CD28", "TNR9")

# Filter out tcsl250 week 1 for combined plots
combined_data <- all_data[!(experiment == "tcsl250" & week == 1)]

# Define a data structure for all CAR sets
car_sets <- list(
  megf8 = list(
    name = "megf8",
    title_prefix = "MEGF8",
    cars = c(cars_to_compare, "MEGF8"),
    time_course_dimensions = c(12, 15),
    dotplot_dimensions = c(14, 6)
  ),
  megf8_all_variants = list(
    name = "megf8_all_variants",
    title_prefix = "MEGF8 Variants",
    cars = c(cars_to_compare, grep("^MEGF8", names(car_colors), value = TRUE)),
    time_course_dimensions = c(15, 5),
    dotplot_dimensions = c(14, 6),
    facet_ncol = 9
  ),
  drb5 = list(
    name = "drb5",
    title_prefix = "DRB5 Group",
    cars = c(cars_to_compare, "DRB5", "DRB528C"),
    time_course_dimensions = c(12, 15),
    dotplot_dimensions = c(14, 6)
  ),
  drb5_all_variants = list(
    name = "drb5_all_variants",
    title_prefix = "All DRB5/DRB528C Variants",
    cars = c(cars_to_compare, grep("^DRB", names(car_colors), value = TRUE)),
    time_course_dimensions = c(18, 12),
    dotplot_dimensions = c(14, 6)
  ),
  tnr9_variants = list(
    name = "tnr9_variants",
    title_prefix = "TNR9 Variants",
    cars = c(cars_to_compare, grep("^TNR9", names(car_colors), value = TRUE)),
    time_course_dimensions = c(18, 12),
    dotplot_dimensions = c(14, 6)
  ),
  cd28_variants = list(
    name = "cd28_variants",
    title_prefix = "CD28 Variants",
    cars = c(cars_to_compare, grep("^CD28", names(car_colors), value = TRUE)),
    time_course_dimensions = c(12, 15),
    dotplot_dimensions = c(14, 6)
  ),
  new_variants = list(
    name = "new_variants",
    title_prefix = "New AI Variants",
    cars = c(cars_to_compare, grep("^COMB", names(car_colors), value = TRUE)),
    time_course_dimensions = c(12, 15),
    dotplot_dimensions = c(14, 6)
  )
)

# Create experiment-specific plots
for (plot_name in names(experiment_plots)) {
  plot_info <- experiment_plots[[plot_name]]
  
  # Filter data based on the expression
  filtered_data <- all_data[eval(plot_info$data_filter)]
  
  # Create the plot
  plots[[plot_name]] <- create_time_course_plot(
    filtered_data,
    title = plot_info$title,
    facet_ncol = plot_info$facet_ncol
  )
}

# Create combined TCSL248 plot
if (!is.null(plots[["tcsl248_wk0"]]) && !is.null(plots[["tcsl248_wk1"]])) {
  plots[["tcsl248_combined"]] <- plots[["tcsl248_wk0"]] / plots[["tcsl248_wk1"]] +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "TCSL248 Incucyte Killing",
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )
  
  # Add dimensions for the combined plot
  experiment_plots[["tcsl248_combined"]] <- list(dimensions = c(8, 20))
}

# Prepare data for dotplots - combine week 0_purified with week 0 and rename experiment
dotplot_data <- copy(combined_data)
# Change tcsl250 with week 0_purified to tcsl250_purified with week 0
dotplot_data[experiment == "tcsl250" & week == "0_purified", `:=`(
  experiment = "tcsl250_purified",
  week = "0"
)]

# Loop through each CAR set and create the plots
for (set_name in names(car_sets)) {
  set_info <- car_sets[[set_name]]

  print("--------------------------------")
  print(set_info$name)
  print(set_info$title_prefix)
  print(unique(set_info$cars))

  
  # Create time course plot
  time_course_key <- paste0("combined_", set_info$name)
  plots[[time_course_key]] <- create_time_course_plot(
    combined_data,
    title = paste("Combined Experiments -", set_info$title_prefix),
    car_subset = unique(set_info$cars),
    #use set_info$ncol if available, else use 4:
    facet_ncol = ifelse(is.null(set_info$facet_ncol), 4, set_info$facet_ncol)
  )
  
  # Create dot plot
  dotplot_key <- paste0("dotplot_", set_info$name)
  plots[[dotplot_key]] <- create_final_timepoint_dotplot(
    dotplot_data,
    title = paste("Final Timepoint Comparison -", set_info$title_prefix),
    car_subset = unique(set_info$cars)
  )
  
  # Store dimensions in experiment_plots
  experiment_plots[[time_course_key]] <- list(dimensions = set_info$time_course_dimensions)
  experiment_plots[[dotplot_key]] <- list(dimensions = set_info$dotplot_dimensions)
}

# Create output directory if it doesn't exist
plots_dir <- file.path(figure_dir, "plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# Create per-plate plots directory
per_plate_dir <- file.path(plots_dir, "per_plate")
dir.create(per_plate_dir, showWarnings = FALSE, recursive = TRUE)

# Function to create a single per-plate diagnostic plot
create_single_plate_plot <- function(data, exp, wk, plt, value_type) {
  # Determine y-axis variable based on value_type
  y_var <- if (value_type == "value_norm") "value_norm" else "value_none"
  y_lab <- if (value_type == "value_norm") "Normalized Intensity to T=0" else "Fraction Cells Remaining"
  
  # Create plot title
  plot_title <- paste0(
    exp, " - ", 
    ifelse(wk == "0_purified", "Week 0 Purified", paste("Week", wk)), 
    " - Plate: ", plt
  )
  
  # Create plot key for the plots list
  plot_key <- paste0(
    "per_plate/", 
    tolower(exp), "_", 
    gsub("_", "", wk), "_", 
    plt, "_", 
    value_type
  )
  
  # Extract unique well-car-ratio combinations
  well_info <- unique(data[, .(well, car_name, ratio, row, col)])
  
  # Create annotation data table
  annotation_dt <- well_info[, .(
    row = row,
    col = col,
    label = paste("Well:", well, 
                 "CAR:", car_name, 
                 "Ratio:", ratio, 
                 sep = "\n"),
    x = min(data$Elapsed) + 2,
    y = 0.9 * max(data[well == well]$value, na.rm = TRUE),
    color = car_colors[car_name]
  )]
    
  # Create the plot
  p <- ggplot(data) +
    geom_line(aes(x = Elapsed, y = value, color = car_name, group=interaction(well, image))) +
    facet_grid(row ~ col, switch = 'y') +
    scale_color_manual(values = car_colors) +
    # Add annotations directly using the data table
    geom_text(data = annotation_dt, 
              aes(x = x, y = y, label = label, color = color),
              hjust = 0, size = 2.5, show.legend = FALSE) +
    labs(
      title = plot_title,
      subtitle = paste0("Metric: ", y_var),
      x = "Elapsed Time (hours)",
      y = y_lab,
      color = "CAR"
    ) +
    theme_minimal() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
      legend.position = "none",  # Remove the legend
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  return(list(
    key = plot_key,
    plot = p,
    dimensions = c(15, 10)
  ))
}

# Get unique experiment-week-plate combinations
plate_combinations <- unique(all_data[!is.na(plate), .(experiment, week, plate)])

# Loop through each plate and value type to create diagnostic plots
for (value_type in c("value_norm", "value_none")) {
  for (i in 1:nrow(plate_combinations)) {
    exp <- plate_combinations[i, experiment]
    wk <- plate_combinations[i, week]
    plt <- plate_combinations[i, plate]
    
    # Filter data for this combination
    plate_data <- all_data[experiment == exp & week == wk & plate == plt & 
                     measurement == "incucyte_count"]
    
    # Skip if no data
    if (nrow(plate_data) == 0) next
    
    # Create the plot
    plot_info <- create_single_plate_plot(plate_data, exp, wk, plt, value_type)
    
    # Add to plots list
    plots[[plot_info$key]] <- plot_info$plot
    experiment_plots[[plot_info$key]] <- list(dimensions = plot_info$dimensions)
    
    cat("Added per-plate plot:", plot_info$key, "\n")
  }
}

# Save all plots
for (name in names(plots)) {
  print(name)
  
  # Handle per-plate plots with subdirectory in the name
  if (grepl("^per_plate/", name)) {
    # Extract the actual filename without the subdirectory prefix
    filename <- gsub("^per_plate/", "", name)
    path <- file.path(plots_dir, "per_plate", paste0(filename, ".pdf"))
  } else {
    path <- file.path(plots_dir, paste0(name, ".pdf"))
  }
  
  print(path)
  
  if (!is.null(plots[[name]]) && name %in% names(experiment_plots)) {
    ggsave(
      filename = path,
      plot = plots[[name]],
      width = experiment_plots[[name]]$dimensions[1],
      height = experiment_plots[[name]]$dimensions[2]
    )
  }
}
