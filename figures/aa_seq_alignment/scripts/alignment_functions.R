#' Amino acid sequence alignment and visualization functions
#' 
#' This file contains functions for aligning amino acid sequences and visualizing
#' the alignments with highlighted differences.

library(Biostrings)
library(stringdist)
library(ggplot2)
library(data.table)
library(grid)
library(gridExtra)
library(msa)

#' Calculate Levenshtein distance between two sequences
#'
#' @param seq1 First sequence
#' @param seq2 Second sequence
#' @return Levenshtein distance
calculate_levenshtein <- function(seq1, seq2) {
  return(stringdist(seq1, seq2, method = "lv"))
}

#' Perform multiple sequence alignment without forcing a specific reference
#'
#' @param sequences Named list of sequences to align
#' @param ref_name Name of the reference sequence (for distance calculation)
#' @return A list containing the aligned sequences and alignment information
align_sequences <- function(sequences, ref_name) {
  # Convert to AAStringSet
  aa_strings <- lapply(sequences, function(seq) AAString(seq))
  aa_set <- AAStringSet(aa_strings)
  names(aa_set) <- names(sequences)
  
  # Perform multiple sequence alignment
  alignment <- msa(aa_set, method = "ClustalW")
  
  # Convert to character strings
  aligned_seqs <- as.character(as(alignment, "AAStringSet"))
  
  # Calculate Levenshtein distances to reference
  ref_seq <- sequences[[ref_name]]
  lev_dists <- sapply(sequences, function(seq) {
    calculate_levenshtein(ref_seq, seq)
  })
  
  return(list(
    aligned_seqs = aligned_seqs,
    lev_dists = lev_dists
  ))
}

#' Create a data frame for alignment visualization
#'
#' @param car_data Data table with CAR sequences
#' @param ref_car Name of the reference CAR
#' @param car_set List of CARs to include in the alignment
#' @return A data frame with alignment information for visualization
prepare_alignment_data <- function(car_data, ref_car, car_set) {
  # Get reference sequence
  ref_seq <- car_data[car_name == ref_car, aa_seq]
  if (length(ref_seq) == 0) {
    stop(paste("Reference CAR", ref_car, "not found in data"))
  }
  
  # Filter car data for the specified car set
  if (!is.null(car_set)) {
    car_names <- car_set$cars
    if (is.null(car_names)) {
      # If NULL, use all CARs
      filtered_data <- car_data
    } else {
      # Filter to include only CARs in the set
      filtered_data <- car_data[car_name %in% car_names]
    }
  } else {
    filtered_data <- car_data
  }
  
  # Ensure reference CAR is included
  if (!(ref_car %in% filtered_data$car_name)) {
    ref_row <- car_data[car_name == ref_car]
    filtered_data <- rbind(ref_row, filtered_data)
  }
  
  # Create named list of sequences for alignment
  sequences <- setNames(
    as.list(filtered_data$aa_seq),
    filtered_data$car_name
  )
  
  # Perform alignment
  alignment_result <- align_sequences(sequences, ref_car)
  
  # Create alignment data table
  alignment_data <- data.table(
    car_name = names(alignment_result$aligned_seqs),
    aligned_seq = unname(alignment_result$aligned_seqs),
    lev_dist = unname(alignment_result$lev_dists)
  )
  
  # Sort by Levenshtein distance
  setorder(alignment_data, lev_dist)
  
  # Add y position (reverse order for plotting)
  alignment_data[, y_pos := .N:1]
  
  return(alignment_data)
}

#' Create a data frame with character-by-character alignment for plotting
#'
#' @param alignment_data Alignment data from prepare_alignment_data
#' @return A data frame with one row per character for plotting
prepare_plot_data <- function(alignment_data) {
  # Get reference sequence and name
  ref_name <- alignment_data[y_pos == max(y_pos), car_name]
  ref_seq <- alignment_data[car_name == ref_name, aligned_seq]
  
  # Split reference sequence into characters
  ref_chars <- strsplit(ref_seq, "")[[1]]
  
  # Create plot data
  plot_data <- alignment_data[, {
    # Split sequence into characters
    chars <- strsplit(aligned_seq, "")[[1]]
    
    # Create data for this sequence
    data.table(
      x_pos = 1:length(chars),
      char = chars
    )
  }, by = .(car_name, y_pos, lev_dist)]
  
  # Add reference character for comparison
  plot_data[, ref_char := ref_chars[x_pos]]
  
  # Add difference flag
  plot_data[, is_diff := char != ref_char]
  
  return(plot_data)
}

#' Create an alignment plot for amino acid sequences
#'
#' @param car_set CAR set from car_sets.R
#' @param ref_car Name of the reference CAR
#' @param car_data Optional data table with CAR sequences (if NULL, will be loaded)
#' @param highlight_color Color for highlighting differences (default: "grey70")
#' @param highlight_alpha Alpha value for highlighting (default: 0.5)
#' @return A ggplot2 object with the alignment visualization
aa_seq_alignment_plot <- function(car_set, ref_car, car_data = NULL, 
                                 highlight_color = "grey90", highlight_alpha = 0.5) {
  # Load CAR data if not provided
  if (is.null(car_data)) {
    car_data <- load_car_data()
  }
  
  # Prepare alignment data
  alignment_data <- prepare_alignment_data(car_data, ref_car, car_set)
  
  # Filter out Zeta if it's just dashes
  if ("Zeta" %in% alignment_data$car_name) {
    zeta_seq <- alignment_data[car_name == "Zeta", aligned_seq]
    if (all(strsplit(zeta_seq, "")[[1]] %in% c("-"))) {
      alignment_data <- alignment_data[car_name != "Zeta"]
      # Recalculate y positions
      alignment_data[, y_pos := .N:1]
    }
  }
  
  # Prepare plot data
  plot_data <- prepare_plot_data(alignment_data)
  
  # Get sequence length and number of sequences for aspect ratio
  seq_length <- max(plot_data$x_pos)
  num_sequences <- length(unique(plot_data$car_name))
  
  # Calculate aspect ratio - make the plot wider than tall
  # This will bring the sequences closer together vertically
  aspect_ratio <- seq_length / (num_sequences * 2)
  
  # Create y-axis labels with row numbers and CAR names
  y_labels <- unique(plot_data[, .(y_pos, car_name)])
  y_labels[, label := sprintf("%d) %s", as.integer(y_pos) - 1, car_name)]
  setorder(y_labels, -y_pos)
  
  # Create plot
  p <- ggplot(plot_data, aes(x = x_pos, y = y_pos)) +
    # Add highlighting for differences
    geom_rect(data = plot_data[is_diff == TRUE], 
              aes(xmin = x_pos - 0.5, xmax = x_pos + 0.5, 
                  ymin = y_pos - 0.3, ymax = y_pos + 0.3),
              fill = highlight_color, alpha = highlight_alpha) +
    # Add characters
    geom_text(aes(label = char, color = car_name), 
              family = "mono", size = 5, vjust = 0.5, fontface = "bold") +
    # Apply CAR colors
    scale_color_manual(values = car_colors) +
    # Set y-axis labels
    scale_y_continuous(
      breaks = y_labels$y_pos,
      labels = y_labels$label,
      expand = c(0.05, 0.05)  # Reduce vertical expansion
    ) +
    # Add spacing between characters
    scale_x_continuous(
      breaks = seq(5, seq_length, by = 5),
      labels = seq(5, seq_length, by = 5)
    ) +
    # Set fixed aspect ratio to bring sequences closer together
    coord_fixed(ratio = aspect_ratio) +
    # Add title
    ggtitle(paste("Alignment to", ref_car)) +
    # Adjust theme
    theme_minimal() +
    theme(
      axis.text.y = element_text(
        family = "mono", 
        size = 10, 
        hjust = 1, 
        color = car_colors[y_labels$car_name],
        face = "bold"
      ),
      axis.text.x = element_text(family = "mono", size = 10),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.margin = margin(20, 20, 20, 20)
    )
  
  return(p)
} 