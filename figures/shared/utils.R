#' Utility functions for loading data and configuration
#' 
#' This file contains utility functions for loading data and configuration
#' that can be used across different scripts.

library(yaml)
library(data.table)

#' Load configuration from YAML file
#'
#' @param config_path Path to the YAML configuration file
#' @return A list with configuration settings
load_config <- function(config_path = "config/paths.yml") {
  if (!file.exists(config_path)) {
    stop(paste("Configuration file not found:", config_path))
  }
  
  config <- yaml::read_yaml(config_path)
  return(config)
}

#' Get path for a data file
#'
#' @param file_name Name of the file
#' @param data_type Type of data (raw, processed, meta)
#' @param config Configuration list from load_config
#' @return Full path to the data file
get_data_path <- function(file_name, data_type = "meta", config = NULL) {
  if (is.null(config)) {
    config <- load_config()
  }
  
  # Determine base path based on data type
  if (data_type == "raw") {
    base_path <- config$cache$raw_data_links
  } else if (data_type == "processed") {
    base_path <- config$cache$processed_data
  } else if (data_type == "meta") {
    base_path <- "data/meta"  # Hardcoded for now, could be added to config
  } else {
    stop(paste("Unknown data type:", data_type))
  }
  
  # Construct full path
  full_path <- file.path(base_path, file_name)
  
  return(full_path)
}

#' Load data from CSV file
#'
#' @param file_name Name of the CSV file
#' @param data_type Type of data (raw, processed, meta)
#' @param config Configuration list from load_config
#' @return A data.table with the loaded data
load_data <- function(file_name, data_type = "meta", config = NULL) {
  full_path <- get_data_path(file_name, data_type, config)
  
  if (!file.exists(full_path)) {
    stop(paste("Data file not found:", full_path))
  }
  
  data <- fread(full_path)
  return(data)
}

#' Load CAR sequence data
#'
#' @param file_name Name of the CSV file with CAR sequences
#' @param config Configuration list from load_config
#' @return A data.table with CAR sequences
load_car_data <- function(file_name = "round2_arrayed_cars.csv", config = NULL) {
  car_data <- load_data(file_name, "meta", config)
  return(car_data)
} 