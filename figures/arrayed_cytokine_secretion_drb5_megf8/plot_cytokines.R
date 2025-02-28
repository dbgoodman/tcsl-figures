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
df_241_242 <- fread(here("data", "processed", "cytokine_secretion", "tcsl241_242_processed.csv"))
df_247 <- fread(here("data", "processed", "cytokine_secretion", "tcsl247_processed.csv"))
df_248 <- fread(here("data", "processed", "cytokine_secretion", "tcsl248_processed.csv"))
df_250 <- fread(here("data", "processed", "cytokine_secretion", "tcsl250_processed.csv"))

# Combine all experiments data for single cytokine plots
df_all_expts <- rbindlist(list(
  df_241_242[, .(expt, car_name, variable = gate_clean, freq)],
  df_247[gate_clean %in% c("IFNg_pos", "IL2_pos", "TNFa_pos"), {
    # Standardize CAR names
    car_name_std = car_name
    car_name_std = ifelse(car_name == "zeta", "Zeta", car_name_std)
    car_name_std = ifelse(car_name == "41BB", "TNR9", car_name_std)
    .(expt, car_name = car_name_std, variable = gate_clean, freq)
  }],
  df_248[gate_clean %in% c("IFNg_pos", "IL2_pos", "TNFa_pos") & subset == "CD3" & CAR == TRUE, 
         .(expt = "TCSL248", car_name, variable = gate_clean, freq)],
  df_250[gate_clean %in% c("IFNg_pos", "IL2_pos", "TNFa_pos") & subset == "CD3" & CAR == TRUE, 
         .(expt = "TCSL250", car_name, variable = gate_clean, freq)]
))

# Filter out KO controls
df_all_expts <- df_all_expts[car_name != "KO"]

# Create combined plot for all experiments
create_combined_plot <- function(df_plot, title_suffix = "") {
  # Define experiment colors
  expt_colors <- c(
    "TCSL241" = "#F8766D",  # red
    "TCSL242" = "#00BA38",  # green
    "TCSL247_D1" = "#619CFF", # blue
    "TCSL247_D2" = "#F564E3", # pink
    "TCSL247_D1V" = "#B79F00", # yellow
    "TCSL248" = "#00BFC4",  # cyan
    "TCSL250" = "#C77CFF"   # purple
  )
  
  # Calculate means for larger points
  df_means <- df_plot[, .(freq = mean(freq)), by=.(expt, variable, car_name)]
  
  # Rank CARs based on means within each experiment and cytokine combination
  df_means[, rank := frank(freq, ties.method = "first"), by=.(expt, variable)]
  df_plot[df_means, rank := i.rank, on=.(expt, variable, car_name)]
  
  # Get number of unique experiments for faceting
  n_expts <- length(unique(df_plot$expt))
  
  # Create plot
  p <- ggplot() +
    # Add small points for replicates
    geom_point(data = df_plot, 
               aes(x = reorder(car_name, rank), y = freq, color = expt),
               size = 1, alpha = 0.5) +
    # Add large points for means
    geom_point(data = df_means,
               aes(x = reorder(car_name, rank), y = freq, color = expt),
               size = 3) +
    facet_wrap(~ expt + variable, 
               scales = "free_y", 
               ncol = 3,
               strip.position = "top") +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.5, "lines"),  # Reduce vertical spacing between facets
      # Only show x-axis text for bottom row
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
    ) +
    scale_color_manual(values = expt_colors) +
    labs(
      x = "CAR",
      y = "% positive cells",
      color = "Experiment",
      title = paste("Individual cytokine secretion by CAR across experiments", title_suffix)
    )
  
  return(p)
}

# Initialize plots list
plots <- list()

# Define gates and CARs
chosen_gate <- c(
  "IL2_pos", "TNFa_pos", "IFNg_pos",
  "IFNg_pos.IL2_pos", "IFNg_pos.TNFa_pos",
  "IL2_pos.TNFa_pos", "IFNg_pos.TNFa_pos.IL2_pos"
)

single_cytokines <- c("IL2_pos", "TNFa_pos", "IFNg_pos")

# Create three different plot versions

# 1. Zeta, CD28, DRB5, TNR9 with TCSL241/242/248/250 (not 247)
zeta_cars <- c("Zeta", "CD28", "DRB5", "TNR9", "MEGF8")
plots[["all_experiments_with_zeta"]] <- create_combined_plot(
  df_all_expts[car_name %in% zeta_cars & 
               expt %in% c("TCSL241", "TCSL242", "TCSL248", "TCSL250")],
  " (All experiments with Zeta)"
)

# 2. CD28, DRB5, DRB528C, TNR9 with all TCSL247 donors and TCSL248/250
drb528c_cars <- c("CD28", "DRB5", "DRB528C", "TNR9")
plots[["tcsl247_248_250_with_drb528c"]] <- create_combined_plot(
  df_all_expts[car_name %in% drb528c_cars & 
               expt %in% c("TCSL247_D1", "TCSL247_D2", "TCSL248", "TCSL250")],
  " (TCSL247/248/250 with DRB528C)"
)

# 3. CD28, DRB5, TNR9 with all experiments and donors
common_cars <- c("CD28", "DRB5", "TNR9")
plots[["all_experiments_common_cars"]] <- create_combined_plot(
  df_all_expts[car_name %in% common_cars & expt != 'TCSL247_D1V'],
  " (All experiments, common CARs)"
)

# 3. All cars in TCSL248/250
hidden_cars <- c("NEG", "POS", "")
plots[["second_round_cars"]] <- create_combined_plot(
  df_all_expts[expt %in% c("TCSL248", "TCSL250") & !(car_name %in% hidden_cars)],
  " (TCSL248/250, all CARs)"
)

# 3. All cytokines in TCSL247/248
hidden_cars <- c("NEG", "POS", "")
plots[["first_round_cars"]] <- create_combined_plot(
  df_all_expts[expt %in% c("TCSL241", "TCSL242") & !(car_name %in% hidden_cars)],
  " (TCSL241/242, all CARs)"
)

# Save plots
plot_dimensions <- list(
  "all_experiments_with_zeta" = c(15, 10),
  "tcsl247_248_250_with_drb528c" = c(15, 10),
  "all_experiments_common_cars" = c(15, 10),
  "second_round_cars" = c(15, 10)

)

# Save each plot
for (name in names(plots)) {
  if (name %in% names(plot_dimensions)) {
    width <- plot_dimensions[[name]][1]
    height <- plot_dimensions[[name]][2]
    
    fn <- here("figures", "arrayed_cytokine_secretion_drb5_megf8",
                     paste0(name, ".pdf"))

    ggsave(
      filename = fn,
      plot = plots[[name]],
      width = width,
      height = height
    )

    message(paste("Saved:\t", fn))
  }

} 
