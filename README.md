
<!-- README.md is generated from README.Rmd. Please edit that file -->

# influenceR

<!-- badges: start -->
<!-- badges: end -->

The `influenceR` function, inside the `influenceR` package, calculates
and evaluates key influence measures, including Cook’s Distance, DFFITS,
and Hadi’s Influence, for linear or generalized linear models by
computing diagnostic metrics to identify influential observations. The
function has an optional setting that will generate plots to visualize
these influence measures, aiding in the detection of data points that
may impact the model’s fit. This tool is helpful for model assessment
and ensuring the reliability of statistical inferences.

## Installation

You can install the development version of influenceR from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("adamglux/A2DATA501", build_vignettes = TRUE)
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(influenceR)
## basic example code

#import dataset mtcars
data1 <- mtcars
model1 <- lm(mpg ~ wt + hp, data = data1)

## run function to create influence measures + plots
sample1 <- influenceR(model1, data1, "all")

## run function to create dataframe only (no plots)
sample2 <- influenceR(model1, data1)

## show dataframe of measures
#sample1$influence_measures

## view / manipulate individual influence measure
#sample2$DFFITs
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" /><img src="man/figures/README-unnamed-chunk-2-2.png" width="100%" /><img src="man/figures/README-unnamed-chunk-2-3.png" width="100%" />
