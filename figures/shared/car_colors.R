#' Shared CAR color definitions
#' 
#' This file contains color definitions for CAR variants that can be used
#' across different figures for consistent visualization.

library(data.table)

# Source utility functions
source("figures/shared/utils.R")

# Load CAR data with colors
car_data <- load_car_data()

# Create a named vector of colors from the car_data
car_colors <- car_data[, setNames(color, car_name)]

# Add colors for controls (not in the CSV)
control_colors <- c(
  "U" = "#999999",
  "K" = "#CCCCCC"
)

# Combine CAR colors with control colors
car_colors <- c(car_colors, control_colors)

#' Get CAR colors
#'
#' @return A named vector of colors for CAR variants
get_car_colors <- function() {
  return(car_colors)
}

#' Apply CAR colors to a ggplot object
#'
#' @param p A ggplot object
#' @param color_aes The aesthetic to apply colors to (color, fill, or both)
#' @return A ggplot object with CAR colors applied
apply_car_colors <- function(p, color_aes = c("color", "fill", "both")) {
  color_aes <- match.arg(color_aes)
  
  if (color_aes == "color" || color_aes == "both") {
    p <- p + scale_color_manual(values = car_colors)
  }
  
  if (color_aes == "fill" || color_aes == "both") {
    p <- p + scale_fill_manual(values = car_colors)
  }
  
  return(p)
} 