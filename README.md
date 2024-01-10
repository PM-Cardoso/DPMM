# DPMM

### Overview

**DPMM** is an R package providing a library of model fitting functions, diagnostics tools for the fitted model and plotting functions. This package is developed from the work done in Cardoso, Dennis, Bowden, Shields and McKinley (2024) <[https://doi.org/10.1186/s12911-023-02400-3](https://doi.org/10.1186/s12911-023-02400-3)>.

### Getting Started

If you are just getting started with **DPMM**, we recommend starting with the tutorial [vignettes](https://pm-cardoso.github.io/DPMM/articles/DPMM.html), the examples throughout the package [documentation](https://pm-cardoso.github.io/DPMM/articles/Worked_Examples.html), and the paper *Dirichlet process mixture models to estimate outcomes for individuals with missing predictor data: application to predict optimal type 2 diabetes therapy in electronic health record data*:

-   Pedro Cardoso, John M. Dennis, Jack Bowden, Beverley M. Shields, Trevelyan J. McKinley the MASTERMIND Consortium. Dirichlet process mixture models to impute missing predictor data in counterfactual prediction models: an application to predict optimal type 2 diabetes therapy. *BMC Medical Informatics and Decision Making* 24, 12 (2024), doi: <[https://doi.org/10.1186/s12911-023-02400-3](https://doi.org/10.1186/s12911-023-02400-3)>.

### Installation

-   Install latest development version from GitHub (requires [devtools](https://github.com/hadley/devtools) package):

``` r
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("PM-Cardoso/DPMM", dependencies = TRUE, build_vignettes = FALSE)
```

This installation won't include the vignettes (they take some time to build), but all of the vignettes are available online at <https://pm-cardoso.github.io/DPMM/index.html>.
