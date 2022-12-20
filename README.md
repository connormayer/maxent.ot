
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

If you publish work that uses `maxent.ot`, please cite this repository.

> `Mayer, C., Tan, A., & Zuraw, K.(2022). maxent.ot: A package for doing Maximum Entropy Optimality Theory in R (Version 0.1.0) [Computer software]. https://doi.org/10.5281/zenodo.7246367`

## Installation

<!--You can install the released version of maxent.ot from [CRAN](https://CRAN.R-project.org) with:

``` r
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
   "extdata", "sample_data_file_small.txt", package = "maxent.ot"
)
# This file has three constraints
data_file_complex <- system.file(
   "extdata", "sample_data_file_large.txt", package = "maxent.ot"
)


# Fit weights to both data sets with simple regularization
simple_model <- optimize_weights(data_file_simple, mu_scalar=0, sigma_scalar=10)
complex_model <- optimize_weights(data_file_complex, mu_scalar=0, sigma_scalar=10)


# Examine predicted probabilities of each model
# Also displayed: log likelihood (of weights given prediction data)
predict_probabilities(data_file_simple, simple_model$weights)
#> $loglik
#> [1] -2.079442
#> 
#> $predictions
#>        UR        SR Freq Constraint1 Constraint3 Predicted Probability
#> 1: Input1 Output1-1    1           1           1                   0.5
#> 2: Input1 Output1-2    1           0           0                   0.5
#> 3: Input2 Output2-1    1           0           1                   0.5
#> 4: Input2 Output2-2    0           0           0                   0.5
#>    Observed Probability Error
#> 1:                  0.5   0.0
#> 2:                  0.5   0.0
#> 3:                  1.0  -0.5
#> 4:                  0.0   0.5
predict_probabilities(data_file_complex, complex_model$weights)
#> $loglik
#> [1] -1.444644
#> 
#> $predictions
#>        UR        SR Freq Constraint1 Constraint2 Constraint3
#> 1: Input1 Output1-1    1           1           0           1
#> 2: Input1 Output1-2    1           0           1           0
#> 3: Input2 Output2-1    1           0           0           1
#> 4: Input2 Output2-2    0           0           1           0
#>    Predicted Probability Observed Probability       Error
#> 1:            0.51385021                  0.5  0.01385021
#> 2:            0.48614979                  0.5 -0.01385021
#> 3:            0.94404417                  1.0 -0.05595583
#> 4:            0.05595583                  0.0  0.05595583


# Compare model fit to training data using the likelihood ratio test
compare_models(simple_model, complex_model, method='lrt')
#>                                     description   chi_sq k_delta   p_value
#> 1 sample_data_file_large~sample_data_file_small 1.269594       1 0.2598428
```

<!--What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so:


```r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/master/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN. -->
