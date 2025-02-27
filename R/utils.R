#' Load configuration from config/paths.yml
#'
#' @return A list containing the configuration
#' @export
load_config <- function() {
  yaml::read_yaml(here::here("config", "paths.yml"))
}

#' Sync data from S3
#'
#' @param direction Either "download" or "upload"
#' @param path Local path relative to data directory
#' @param config Optional config list (will be loaded if not provided)
#' @export
sync_s3 <- function(direction = c("download", "upload"), path = NULL, config = NULL) {
  direction <- match.arg(direction)
  if (is.null(config)) config <- load_config()
  
  # Build S3 URI
  bucket <- config$data_sources$s3$bucket
  prefix <- config$data_sources$s3$prefix
  s3_uri <- sprintf("s3://%s/%s", bucket, prefix)
  
  # Build local path
  local_path <- if (is.null(path)) {
    here::here("data")
  } else {
    here::here("data", path)
  }
  
  # Build aws s3 sync command
  cmd <- if (direction == "download") {
    sprintf("aws s3 sync %s %s", s3_uri, local_path)
  } else {
    sprintf("aws s3 sync %s %s", local_path, s3_uri)
  }
  
  # Execute command
  system(cmd)
}

#' Get path to analysis repo
#'
#' @param ... Additional path components to append
#' @param config Optional config list (will be loaded if not provided)
#' @return Character string with full path
#' @export
analysis_path <- function(..., config = NULL) {
  if (is.null(config)) config <- load_config()
  file.path(config$data_sources$analysis_repo$root, ...)
}

#' Get path to Box data
#'
#' @param ... Additional path components to append
#' @param config Optional config list (will be loaded if not provided)
#' @return Character string with full path
#' @export
box_path <- function(..., config = NULL) {
  if (is.null(config)) config <- load_config()
  file.path(config$data_sources$box$root, ...)
} 