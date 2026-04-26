# R code accompanying the manuscript *Estimating the risk of sustained transmission from early infectious disease outbreak data*

This repository contains the R scripts used to estimate the probability of a major outbreak (PMO) from early incidence data via three methods, analytic, trajectory matching (TM) and machine learning (ML), as described in the manuscript.

More specifically, the repository supports:

1. simulation of outbreak datasets under branching-process transmission models
2. training of an ML classifier on simulated data to estimate PMO
3. computation of PMO for incidence datasets analytically and via TM/ML
4. computation of PMO under underreporting

## Main scripts

### Core data generation and modelling

- `serial_interval_discretisation.R`  
  Generates discretised serial interval weights for downstream use.

- `simulate_training_dataset.R`  
  Simulates training data under a chosen prior on `I_1` and saves a dataset for TM and ML estimation.

- `compute_analytic_pmo.R`  
  Contains functions for analytic PMO estimation and computes analytic PMO values over a set of incidence datasets.
  
- `training_and_prediction_TM_and_ML.R`  
  Trains an ML classifier on a chosen simulated dataset and predicts PMO for a set of incidence datasets using TM and ML methods.

### Underreporting pipeline

- `merge_outputs_and_likelihood.R`  
  Merges analytic, TM and ML PMO outputs over incidence datasets and computes the dataset probabilities.

- `calculate_pmo_underreporting.R`  
  Applies the binomial reporting model for a chosen observed dataset `Y`, updates dataset probabilities to posterior weights and saves outputs required for plotting.

- `plot_underreporting_pmo.R`  
  Produces underreporting plots, including posterior bar plots over the true underlying incidence datasets and PMO comparison panels across reporting probabilities.

### Convergence and comparison plots

- `training_ML_convergence.R`  
  Trains ML models at increasing simulation sizes and saves outputs for convergence analysis.

- `plot_convergence.R`  
  Plots convergence of ML and TM PMO estimates against the analytic benchmark for selected incidence datasets.

- `plot_pmo_comparison_summary.R`  
  Produces PMO comparison summaries for selected incidence datasets at a fixed training size.

### Other plotting / simulation scripts

- `simplot.R`  
  Generates simulation-based plots used in the manuscript (Figure 1).

## Expected workflow

A typical workflow for the main 3-week analysis is:

1. run `serial_interval_discretisation.R`
2. run `simulate_training_dataset.R`
3. run `training_and_prediction_TM_and_ML.R`
4. run `compute_analytic_pmo.R`
5. run `merge_outputs_and_likelihood.R`
6. run `calculate_pmo_underreporting.R`
7. run `plot_underreporting_pmo.R`

For convergence analyses:

1. run `training_ML_convergence.R`
2. run `plot_convergence.R`
3. run `plot_pmo_comparison_summary.R`

## Key outputs

Main generated outputs:

- discretised serial interval weights
- simulated training datasets
- analytic, TM, and ML PMO predictions
- PMO summary tables across reporting probabilities
- manuscript plots

## Software

The code is written in R and primarily uses packages including:

- `data.table`
- `dplyr`
- `tidyr`
- `ggplot2`
- `ranger`
- `ggnewscale`