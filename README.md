# DPMM

### Overview

**DPMM** is an R package providing a library of model fitting functions, diagnostics tools for the fitted model and plotting functions. This package is developed from the work done in the preprint Cardoso, Dennis, Bowden, Shields and McKinley (2023) <https://doi.org/10.1101/2022.07.26.22278066>.

### Getting Started

If you are just getting started with **DPMM**, we recommend starting with the tutorial [vignettes](https://www.google.com), the examples throughout the package [document](https://www.google.com), and the paper *Dirichlet process mixture models to estimate outcomes for individuals with missing predictor data: application to predict optimal type 2 diabetes therapy in electronic health record data*:

-   Pedro Cardoso, John M. Dennis, Jack Bowden, Beverley M. Shields, Trevelyan J. McKinley the MASTERMIND Consortium. Dirichlet process mixture models to estimate outcomes for individuals with missing predictor data: application to predict optimal type 2 diabetes therapy in electronic health record data. *medRxiv*, doi: <https://doi.org/10.1101/2022.07.26.22278066>. ([medRxiv](https://doi.org/10.1101/2022.07.26.22278066))

### Installation

-   Install latest development version from GitHub (requires [devtools](https://github.com/hadley/devtools) package):

``` r
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("PM-Cardoso/DPMM", dependencies = TRUE, build_vignettes = FALSE)
```

This installation won't include the vignettes (they take some time to build), but all of the vignettes are available online [here](https://www.google.com).
