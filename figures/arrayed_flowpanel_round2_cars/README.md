# Arrayed Flow Panel Round 2 CAR Analysis

This directory contains scripts for processing and analyzing flow cytometry data from the TCSL248 and TCSL250 experiments.

## Directory Structure

- `main.R`: Main script that processes both experiments and generates individual and combined plots
- `process_experiment.R`: Contains a function to process a single experiment
- `data_setup/`: Contains core data processing and plotting functions
  - `process_flowpanel_data.R`: Functions for processing flow cytometry data
  - `plot_flowpanel.R`: Functions for generating plots from processed data

## Usage

To process both experiments and generate all plots, run:

```bash
Rscript figures/arrayed_flowpanel_round2_cars/main.R
```

## Output

The script generates the following outputs:

1. Processed data files:
   - `data/processed/flowpanel_round2/tcsl248_processed.rds`
   - `data/processed/flowpanel_round2/tcsl250_processed.rds`

2. Individual experiment plots:
   - `figures/arrayed_flowpanel_round2_cars/plots/tcsl248/`
   - `figures/arrayed_flowpanel_round2_cars/plots/tcsl250/`

3. Combined experiment plots:
   - `figures/arrayed_flowpanel_round2_cars/plots/combined/`

## Experiment Details

### TCSL248
- Contains unstimulated (X) and stimulated (A, B, C) replicates
- Stimulation ratios: 1:8 (A, B) and 1:4 (C)

### TCSL250
- Contains only stimulated replicates (A, B, C)
- All replicates use 1:4 stimulation ratio 