
<!-- README.md is generated from README.Rmd. Please edit that file -->

# maxent.ot

<!-- badges: start -->

[![DOI](https://zenodo.org/badge/272315147.svg)](https://zenodo.org/badge/latestdoi/272315147)
<!-- badges: end -->

This package allows you to fit Maximum Entropy Optimality Theory models
to data sets, generate the predictions made by such models for novel
data, and compare the fit of different models using a variety of
metrics. This package is still in development, and is being prepared for
submission to CRAN.

The authors of this package are [Connor Mayer](http://connormayer.com),
[Adeline Tan](https://linguistics.ucla.edu/person/adeline-tan/), and
[Kie Zuraw](https://linguistics.ucla.edu/people/zuraw/).

## Citing maxent.ot

If you publish work that uses `maxent.ot`, please cite the following
paper and repository:

> `Mayer, C., Tan, A., & Zuraw, K. (in press). Introducing maxent.ot: an R package for Maximum Entropy constraint grammars. Phonological Data and Analysis.`

> `Mayer, C., Tan, A., & Zuraw, K.(2024). maxent.ot: A package for doing Maximum Entropy Optimality Theory in R (Version 1.0.0) [Computer software]. 10.5281/zenodo.7246366`

## Installation

<!--You can install the released version of maxent.ot from [CRAN](https://CRAN.R-project.org) with:
  &#10;  ``` r
install.packages("maxent.ot")
```
-->

You can install the development version from
[GitHub](https://github.com/) with:

``` r
if (!require(devtools)) {
install.packages("devtools", repos = "http://cran.us.r-project.org")
}
if (!require(maxent.ot)) {
devtools::install_github("connormayer/maxent.ot")
}
```

## Example

This is a simple example workflow of fitting two MaxEnt OT models to the
same data (with different constraint sets), examining their predicted
frequencies, and comparing their fits using the likelihood ratio test.

``` r
library(maxent.ot)

# Get paths to input files.
# This file has two constraints
data_file_simple <- system.file(
  "extdata", "sample_data_frame.csv", package = "maxent.ot"
)
# This file has three constraints
data_file_complex <- system.file(
  "extdata", "sample_data_frame_large.csv", package = "maxent.ot"
)

# Read files into data frames
df_simple <- read.csv(data_file_simple)
df_complex <- read.csv(data_file_complex)

# Fit weights to both data sets with simple regularization
simple_model <- optimize_weights(df_simple, mu=0, sigma=10)
complex_model <- optimize_weights(df_complex, mu=0, sigma=10)

# Examine predicted probabilities of each model
# Also displayed: log likelihood (of weights given prediction data)
predict_probabilities(df_simple, simple_model$weights)
#> $loglik
#> [1] -1.444645
#> 
#> $predictions
#>    Input    Output Freq Constraint1 Constraint2  Predicted Observed       Error
#> 1 Input1 Output1-1    1           1           0 0.51384754      0.5  0.01384754
#> 2 Input1 Output1-2    1           0           1 0.48615246      0.5 -0.01384754
#> 3 Input2 Output2-1    1           0           0 0.94404279      1.0 -0.05595721
#> 4 Input2 Output2-2    0           0           1 0.05595721      0.0  0.05595721
predict_probabilities(df_complex, complex_model$weights)
#> $loglik
#> [1] -1.444644
#> 
#> $predictions
#>    Input    Output Freq Constraint1 Constraint2 Constraint3  Predicted Observed
#> 1 Input1 Output1-1    1           1           0           1 0.51385019      0.5
#> 2 Input1 Output1-2    1           0           1           0 0.48614981      0.5
#> 3 Input2 Output2-1    1           0           0           1 0.94404422      1.0
#> 4 Input2 Output2-2    0           0           1           0 0.05595578      0.0
#>         Error
#> 1  0.01385019
#> 2 -0.01385019
#> 3 -0.05595578
#> 4  0.05595578


# Compare model fit to training data using the likelihood ratio test
compare_models(simple_model, complex_model, method='lrt')
#>            description       chi_sq k_delta   p_value
#> 1 df_complex~df_simple 2.451046e-06       1 0.9987508
```
