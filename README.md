# TCSL Figures Repository

This repository contains the code for generating figures and final analysis for the TCSL paper. It works in conjunction with the main analysis repository at `/Users/dbg/code/tcsl`.

## Repository Structure

- `R/`: R functions and utilities
- `config/`: Configuration files including data paths
- `data/`: Data directory (see Data Management below)
- `output/`: Generated figures and tables
- `tests/`: Test cases for R functions
- `figures/`: R scripts for individual figures

## Setup and Dependencies

```r
# Install renv for dependency management
install.packages("renv")
renv::init()

# Install required packages
renv::install()
```

## Data Management

This repository works with data from multiple sources:
1. Main analysis repository: `/Users/dbg/code/tcsl`
2. Box filesystem: `/Users/dbg/Library/CloudStorage/Box-Box/tcsl`
3. S3 storage for large datasets

### Data Organization
- Large datasets are stored in S3
- Box data is accessed directly from the Box sync directory
- Processed data for figures is cached in `data/processed/`

### Working with S3 Data

The repository includes utility functions for syncing with S3:

```r
# Download data from S3
sync_s3("download")

# Upload data to S3
sync_s3("upload")
```

Make sure you have AWS CLI configured with appropriate credentials.

## Figure Generation

Each figure has its own R script in the `figures/` directory. To generate a specific figure:

```r
source("figures/figure1.R")
```

## Helper Functions

The `R/utils.R` file contains helper functions for working with paths:

```r
# Get path to analysis repo
analysis_path("some", "path")

# Get path to Box data
box_path("some", "path")
```

## Version Information

- Main Analysis Repository: `/Users/dbg/code/tcsl`
- Box Data Path: `/Users/dbg/Library/CloudStorage/Box-Box/tcsl`
- S3 Bucket: `tcsl-data`

## Contributing

[Add contribution guidelines]

## License

[Add license information] 