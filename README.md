# DPMM_fit

###
Fitting a DPMM to complete/incomplete datasets with continuous and categorical variables.

---
Datasets included in this repository:
  - Dataset 1: two continuous variables with distinct clusters of data
  - Dataset 2: two continuous variables + one categorical variable. Two separate clusters (same X1 values), only distinction between both clusters is the categorical variable
  - Dataset 3: two categorical variables with distint clusters of data
  - Dataset 4: Y variables, three continuous + one categorical predictor variables
  - Dataset 5: two continuous variables + two categorical variables. Two separate clusters (same X1 values), only distinction between both clusters is the categorical variables
  - Dataset 6: two continuous variables with distinct clusters of data. Equivalent to Dataset 1, but more correlation between both variables

---

Fitting the model is done with `dpmm_fit.R`. This file includes the function _runModel(dataset, mcmc_iterations, L, standardise)_:
  - dataset: a dataframe with any number of continuous and categorical variables with/without missing values
  - mcmc_iterations: number of iterations of the MCMC algorithm
  - L: maximum number of components to fit the model (should be below optimal number of components)
  - standardise: logical (default = TRUE), standardises all continuous variables
Returns a list of class _dpmm_fit_ with all the information provided and standardisation terms.

`conditional_RW.R` and `conditional_RW_block.R` provide conditional gaussian updates for missing continuous variables.

---
#### Plots

Traceplots for all DPMM parameters can be plotted through `plot_dpmm_fit.R`. This includes the function _plot.dpmm_fit(x, trace = TRUE, density = TRUE)_:
  - x: an object of class _dpmm_fit_
  - trace: logical (default = TRUE) plot traceplots
  - density: logical (default = TRUE) plot density plots
Returns a series of plots for all parameters.

The number of components used in the DPMM can be analysed through `plot_alpha.R`. This compares several DPMMs by plotting the number of components with used through the iterations, the average number of individuals in ranked components and a traceplot of alpha values. This file includes the function _plot.alpha(x)_:
  - x: a list of several _dpmm_fit_ objects
Returns a plot divided into three sections.

A plot of random samples can be generated with `plot_ggpairs.R`. This can produce a plot for DPMM samples or a comparison versus a new dataset (with equal number of samples) when new dataset provided. This file has the function _plot.ggpairs(x, newdata, iterations, nburn)_:
  - x: object of class _dpmm_fit_ or list of _dpmm_fit_ objects
  - newdata: a dataframe with identical structure to data in DPMM fit. _newdata_ can only be used with _nburn_. Random samples from _x_ will have the same number of draws as _newdata_.
  - iterations: vector of iterations for random samples of _x_. This can only be used with _x_ and will provide a GGally::ggpairs() plot for random samples of _x_.
  - nburn: cut-off point for burning iterations. Can only be used when _newdata_ is supplied, and cannot be supplied alongside _iterations_
Returns a generalised pairs plot.

