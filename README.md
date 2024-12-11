# transferlearning


This repository contains the code and data used in the paper:

**"Transfer learning for piecewise-constant mean estimation: Optimality, $\ell_1$- and $\ell_0$-penalisation"**  
Authors: Fan Wang and Yi Yu

## Overview

The repository supports all experiments and analyses presented in the paper. The implementation is done in R (version >= 4.1).

## Dependencies

To reproduce the results, the following R packages are required:

```r
list.of.packages <- c("genlasso", "changepoints", "latex2exp", "egg", "tictoc", "ggplot2", "plotrix")
new_packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
for (package in new_packages) {
  install.packages(package, repos = "http://cran.us.r-project.org", quiet = TRUE)
}
```

## File Descriptions

Below is a summary of the key scripts and their purposes:

| **File Name**                          | **Description**                                                                                     |
|----------------------------------------|-----------------------------------------------------------------------------------------------------|
| `transfer_learning_scenario_1_original.R` | Simulation studies for Scenario 1 (Section 5.1): Equally spaced change points.                     |
| `transfer_learning_scenario_2_original.R` | Simulation studies for Scenario 2 (Section 5.1): Unequally spaced change points.                   |
| `transfer_learning_add.R`              | Simulation studies with varied source frequencies (Section 5.1.1).                                 |
| `transfer_learning_extensions.R`       | Simulation studies incorporating target data (Section 5.1.2).                                      |
| `transfer_learning_all.R`              | Generates Figure 1 in the paper.                                                                   |
| `transfer_learning_tk.R`               | Sensitivity analysis on the screening size.                                                        |
| `transfer_learning_tempde.R`           | Simulations with dependent noise variables and discrepancy dependence between target and sources.   |
| `transfer_learning_large_size.R`       | Simulation studies for varying sample sizes (`n0`).                                                |
| `transfer_learning_GDP.R`              | Real data analysis for the GDP and IP dataset.                                                     |
| `transfer_learning_US_electricity.R`   | Real data analysis for the U.S. electricity dataset.                                               |
| `transfer_learning_air_quality.R`      | Real data analysis for the air quality dataset.                                                    |
| `transfer_learning_functions.R`        | Contains all functions used across methods.                                                        |

## Datasets

The repository also includes several datasets:

- **Air Quality**: Air quality data for different cities.  
- **`gdp_q.csv`**: Quarterly GDP dataset.  
- **`ip.csv`**: Monthly Industrial Production (IP) dataset.  




