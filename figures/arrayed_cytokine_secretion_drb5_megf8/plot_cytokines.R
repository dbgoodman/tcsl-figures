library(data.table)
library(ggplot2)
library(patchwork)
library(stringr)
library(here)

# Function to prepare data for plotting
prepare_data <- function(df_combined, selected_cars = NULL, gates = NULL) {
  # Base filtering
  dt <- df_combined[
    CAR == T & 
    !(car_name %in% c('POS','NEG')) & 
    !is.na(car_name) & 
    car_name != "" & 
    subset == "CD3"
  ]
  
  # Filter by selected cars if provided
  if (!is.null(selected_cars)) {
    dt <- dt[car_name %in% selected_cars]
  }
  
  # Filter by gates if provided
  if (!is.null(gates)) {
    dt <- dt[gate_clean %in% gates]
  }
  
  return(dt)
}

# Function to create heatmap
create_heatmap <- function(data, gates, is_relative = TRUE, title_suffix = "") {
  # Calculate mean frequencies and z-scores
  plot_data <- data[
    gate_clean %in% gates, 
    .(mean_freq = mean(freq)), 
    by=.(car_name, subset, gate_clean, expt)
  ][
    , scaled_mean_freq := scale(mean_freq), 
    by=.(subset, gate_clean, expt)
  ][
    , car_name := reorder(car_name, scaled_mean_freq)
  ][
    , gate_clean := factor(gate_clean, 
                          levels = unique(gate_clean[order(str_count(gate_clean, "\\."))]))
  ]
  
  # Create plot
  p <- ggplot(plot_data) +
    geom_tile(aes(y=car_name, 
                  x=gate_clean, 
                  fill=if(is_relative) scaled_mean_freq else mean_freq)) +
    facet_grid(~expt) +
    scale_fill_gradientn(
      colors = rev(hcl.colors(100, "Peach"))) +
    geom_text(aes(y=car_name, 
                  x=gate_clean, 
                  label=round(mean_freq)), 
              size=3) +
    scale_x_discrete(labels=function(x) gsub('\\.', '\n',gsub('_pos','',x))) +
    labs(x='Cytokine combinations', 
         y='CARs',
         fill=if(is_relative) 
           'Relative cytokine secretion\nz-score xform\n% positive labelled' 
           else 'Absolute cytokine secretion\n% positive',
         title=paste0(if(is_relative) 
           'Relative' else 'Absolute', 
           ' cytokine secretion across TCSL248 and TCSL250', 
           title_suffix))
  
  return(p)
}

# Function to create scatter plot
create_scatter_plot <- function(data, single_cytokines, title_suffix = "", show_replicates = FALSE) {
  # Calculate frequencies and create combined labels
  plot_data <- data[
    gate_clean %in% single_cytokines
  ]
  
  if (show_replicates) {
    # Create mean data for larger points
    mean_data <- plot_data[
      , .(freq = mean(freq)), 
      by=.(car_name, gate_clean, expt)
    ][
      , car_name_subset := paste(car_name, gate_clean, sep=".")
    ]
    
    # Add car_name_subset to replicate data for consistent ordering
    plot_data[, car_name_subset := paste(car_name, gate_clean, sep=".")]
  } else {
    # Just use means as before
    plot_data <- plot_data[
      , .(freq = mean(freq)), 
      by=.(car_name, gate_clean, expt)
    ][
      , car_name_subset := paste(car_name, gate_clean, sep=".")
    ]
  }
  
  # Create base plot
  p <- ggplot() +
    facet_wrap(~gate_clean, scales="free",
               labeller=labeller(gate_clean=function(x) gsub('_pos', '', x))) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_color_brewer(palette="Set1") +
    scale_x_discrete(labels=function(x) gsub('\\.[^.]+$', '', x)) +
    labs(x='CAR, ranked separately per cytokine', 
         y='% positive cells', 
         color='Experiment',
         title=paste0('Individual cytokine secretion by CAR across experiments', 
                     title_suffix))
  
  if (show_replicates) {
    # Add small points for replicates
    p <- p + 
      geom_point(data=plot_data,
                 aes(x=reorder(car_name_subset, freq), y=freq, color=expt),
                 size=1, alpha=0.5) +
      # Add large points for means
      geom_point(data=mean_data,
                 aes(x=reorder(car_name_subset, freq), y=freq, color=expt),
                 size=3)
  } else {
    # Original single point behavior
    p <- p + 
      geom_point(data=plot_data,
                 aes(x=reorder(car_name_subset, freq), y=freq, color=expt),
                 size=3)
  }
  
  return(p)
}

# Load and prepare data
df_248 <- fread(here('data', 'processed', 'cytokine_secretion', 'tcsl248_processed.csv'))
df_250 <- fread(here('data', 'processed', 'cytokine_secretion', 'tcsl250_processed.csv'))

df_248[, expt := 'TCSL248']
df_250[, expt := 'TCSL250']

df_combined <- rbind(df_248, df_250)

# Define gates and CARs
chosen_gate <- c(
  "IL2_pos", "TNFa_pos", "IFNg_pos",
  "IFNg_pos.IL2_pos", "IFNg_pos.TNFa_pos",
  "IL2_pos.TNFa_pos", "IFNg_pos.TNFa_pos.IL2_pos"
)

single_cytokines <- c("IL2_pos", "TNFa_pos", "IFNg_pos")
selected_cars <- c("DRB5", "MEGF8", "CD28", "TNR9", "Zeta")

# Prepare datasets
data_all <- prepare_data(df_combined)
data_filtered <- prepare_data(df_combined, selected_cars)

# Create plots
plots <- list()

# Full dataset plots
plots[["relative_cytokine_secretion"]] <- create_heatmap(data_all, chosen_gate, is_relative = TRUE)
plots[["absolute_cytokine_secretion"]] <- create_heatmap(data_all, chosen_gate, is_relative = FALSE)
plots[["individual_cytokines"]] <- create_scatter_plot(data_all, single_cytokines, show_replicates = TRUE)

# Filtered dataset plots
plots[["relative_cytokine_secretion_filtered"]] <- create_heatmap(data_filtered, chosen_gate, 
                                                                 is_relative = TRUE, " (Selected CARs)")
plots[["absolute_cytokine_secretion_filtered"]] <- create_heatmap(data_filtered, chosen_gate, 
                                                                 is_relative = FALSE, " (Selected CARs)")
plots[["individual_cytokines_filtered"]] <- create_scatter_plot(data_filtered, single_cytokines, 
                                                              " (Selected CARs)", show_replicates = TRUE)

# Save plots
plot_dimensions <- list(
  "relative_cytokine_secretion" = c(15, 7),
  "absolute_cytokine_secretion" = c(15, 7),
  "individual_cytokines" = c(12, 5),
  "relative_cytokine_secretion_filtered" = c(15, 5),
  "absolute_cytokine_secretion_filtered" = c(15, 5),
  "individual_cytokines_filtered" = c(12, 5)
)

# Save all plots
for (name in names(plots)) {
  ggsave(
    here('figures', 'arrayed_cytokine_secretion_drb5_megf8', paste0(name, '.pdf')),
    plots[[name]], 
    width = plot_dimensions[[name]][1],
    height = plot_dimensions[[name]][2]
  )
} 