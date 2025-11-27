# CODE for "Simplicity versus Complexity: The Role of Historical Average in Kelly, Malamud, and Zhou's (2024) RFF Model"

This is the ReadMe document for running the real data analysis and results presented in the paper.

## Summary of the paper

Under the null of zero slope coefficients, Kelly, Malamud, and Zhou’s (2024) random Fourier features (RFF) model, when estimated at the optimal regularization level, converges asymptotically to an intercept-only specification. Consistent with this null, the RFF model with an intercept matches the historical average across sample periods and training-window lengths. This benchmark also explains the authors’ preferred but misspecified no-intercept version, which embeds an implicit intercept from non-centered predictors. When predictors are mean-centered, shutting off this approximation mechanism, the no-intercept RFF model exhibits negligible market-timing ability. Moreover, neither the benchmark nor the RFF model yields significant CAPM $\alpha$ estimates after 1975.
