#' Shared CAR sets definitions
#' 
#' This file contains definitions for CAR sets that can be used
#' across different figures for consistent visualization.

# Source required files
source("figures/shared/utils.R")
source("figures/shared/car_colors.R")

# Define base CARs to compare against
cars_to_compare <- c("Zeta", "CD28", "TNR9")

# Define CAR sets for different analyses
car_sets <- list(
  all_cars = list(
    name = "all_cars",
    title_prefix = "All CARs",
    cars = NULL,  # NULL means include all CARs
    dimensions = c(12, 8)
  ),
  megf8 = list(
    name = "megf8",
    title_prefix = "MEGF8",
    cars = c("MEGF8"),
    dimensions = c(12, 8)
  ),
  megf8_all_variants = list(
    name = "megf8_all_variants",
    title_prefix = "MEGF8 Variants",
    cars = c( "MEGF8", grep('^MEGF8\\.', names(car_colors), value = TRUE)),
    dimensions = c(12, 8)
  ),
  drb5 = list(
    name = "drb5",
    title_prefix = "DRB5 Group",
    cars = c("DRB5", "DRB528C"),
    dimensions = c(12, 8)
  ),
  drb5_all_variants = list(
    name = "drb5_all_variants",
    title_prefix = "All DRB5/DRB528C Variants",
    cars = c("DRB5", "DRB528C", grep('^DRB5.*\\.', names(car_colors), value = TRUE)),
    dimensions = c(12, 8)
  ),
  tnr9_variants = list(
    name = "tnr9_variants",
    title_prefix = "TNR9 Variants",
    cars = c(grep('^TNR9', names(car_colors), value = TRUE)),
    dimensions = c(12, 8)
  ),
  cd28_variants = list(
    name = "cd28_variants",
    title_prefix = "CD28 Variants",
    cars = c(grep('^CD28', names(car_colors), value = TRUE)),
    dimensions = c(12, 8)
  ),
  new_variants = list(
    name = "new_variants",
    title_prefix = "New AI Variants",
    cars = c(grep('^COMB', names(car_colors), value = TRUE)),
    dimensions = c(12, 8)
  )
)

#' Get CAR sets
#'
#' @return A list of CAR sets
get_car_sets <- function(cars_to_compare = c('CD28','TNR9','Zeta')) {
  # Create a deep copy of the car_sets list
  modified_car_sets <- car_sets
  
  # Add the cars_to_compare to each car set
  for (set_name in names(modified_car_sets)) {
    # For the all_cars set, we don't need to modify the cars list
    if (set_name != "all_cars" && !is.null(modified_car_sets[[set_name]]$cars)) {
      # Add the cars_to_compare to the existing cars
      modified_car_sets[[set_name]]$cars <- unique(c(cars_to_compare, modified_car_sets[[set_name]]$cars))
    }
  }
  
  return(modified_car_sets)
}

#' Filter data for a specific CAR set
#'
#' @param data The data to filter
#' @param car_set The CAR set to filter for
#' @param car_column The column containing CAR names (default: "car_name")
#' @return Filtered data containing only CARs in the specified set
filter_car_set <- function(data, car_set, car_column = "car_name") {
  # If cars is NULL, return all data
  if (is.null(car_set$cars)) {
    return(data)
  }
  
  # Filter data to include only CARs in the set
  filtered_data <- data[get(car_column) %in% car_set$cars]
  
  return(filtered_data)
} 