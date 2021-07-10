
<!-- README.md is generated from README.Rmd. Please edit that file -->

# maxent.ot

<!-- badges: start -->
<!-- badges: end -->

This package allows you to fit Maximum Entropy Optimality Theory models
to data sets, generate the predictions made by such models for novel
data, and compare the fit of different models using a variety of
metrics. This package is still in development, and is being prepared for
submission to CRAN.

The authors of this package are [Connor Mayer](http://connormayer.com)
and [Kie Zuraw](https://linguistics.ucla.edu/people/zuraw/).

## Installation

<!--You can install the released version of maxent.ot from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("maxent.ot")
```
-->

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("connormayer/maxent.ot")
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
   "extdata", "sample_data_file.txt", package = "maxent.ot"
)
# This file has three constraints
data_file_complex <- system.file(
   "extdata", "sample_data_file_2.txt", package = "maxent.ot"
)

# Fit weights to both data sets with simple regularization
simple_model <- optimize_weights(data_file_simple, mu_scalar=0, sigma_scalar=10)
complex_model <- optimize_weights(data_file_complex, mu_scalar=0, sigma_scalar=10)

# Examine predicted probabilities of each model
simple_predictions <- predict_probabilities(data_file_simple, simple_model$weights)
complex_predictions <- predict_probabilities(data_file_complex, complex_model$weights)

# Compare model fit to training data using the likelihood ratio test
compare_models(simple_model, complex_model, method='lrt')
#>                           description    chi_sq k_delta   p_value
#> 1 sample_data_file_2~sample_data_file 0.1334752       1 0.2851443
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
