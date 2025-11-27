# CODE for "Simplicity versus Complexity: The Role of Historical Average in Kelly, Malamud, and Zhou's (2024) RFF Model"

This is the ReadMe document for running the real data analysis and results presented in the paper.

## Summary of the paper

Under the null of zero slope coefficients, Kelly, Malamud, and Zhou’s (2024) random Fourier features (RFF) model, when estimated at the optimal regularization level, converges asymptotically to an intercept-only specification. Consistent with this null, the RFF model with an intercept matches the historical average across sample periods and training-window lengths. This benchmark also explains the authors’ preferred but misspecified no-intercept version, which embeds an implicit intercept from non-centered predictors. When predictors are mean-centered, shutting off this approximation mechanism, the no-intercept RFF model exhibits negligible market-timing ability. Moreover, neither the benchmark nor the RFF model yields significant CAPM $\alpha$ estimates after 1975.

## 1. Predictions

The **`Step1_Predictions`** folder contains all code required to generate monthly stock return forecasts.

### Files and Functions

- **`GYdata.mat`**  
  Contains the raw predictor variables and equity returns from the Goyal and Welch (2008) dataset.

- **`rff_main.m`**  
  Implements Random Fourier Feature (RFF) construction and time-series prediction.  
  The function supports:
  - **Fixed λ ridge regression**
  - **Optimally tuned λ ridge regression**

  You may adjust parameter settings to evaluate two scenarios:
  1. `demeanX = 0`, `demeanY = 0` — *no intercept (raw features and returns)*  
  2. `demeanX = 1`, `demeanY = 1` — *with intercept (demeaned features and returns)*

- **`GW_benchmark_main.m`**  
  Implements benchmark prediction models following Goyal & Welch (2008).  
  Similar to `rff_main.m`, the function provides:
  - **Fixed λ ridge regression**
  - **Optimal λ ridge regression**

  Available configurations:
  - `demeanX = 0`, `demeanY = 0` (no intercept)  
  - `demeanX = 1`, `demeanY = 1` (with intercept)

- **`Nagel_2025.m`**  
  Performs performance comparisons consistent with Nagel (2025).

## 2. Exhibits

The **`Step2_Exhibits`** folder contains scripts used to construct tables and summary exhibits based on the prediction outputs.

### Files and Functions

- **`Merge_data.R`**  
  Merges output from different RFF configurations and computes averaged predictive performance.

- **`table.R`**  
  Generates the main prediction performance table for the paper.

- **`table_Nagel.R`**  
  Replicates and reports benchmark results following Nagel (2025).
