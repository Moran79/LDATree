
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LDATree <a href="http://iamwangsiyu.com/LDATree/"><img src="man/figures/logo.png" align="right" height="139" alt="LDATree website" /></a>

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/LDATree)](https://CRAN.R-project.org/package=LDATree)
[![R-CMD-check](https://github.com/Moran79/LDATree/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Moran79/LDATree/actions/workflows/R-CMD-check.yaml)
![CRAN Downloads](https://cranlogs.r-pkg.org/badges/grand-total/LDATree)
<!-- badges: end -->

`LDATree` is an R modeling package for fitting classification trees. If
you are unfamiliar with classification trees, here is a
[tutorial](http://www.sthda.com/english/articles/35-statistical-machine-learning-essentials/141-cart-model-decision-tree-essentials/)
about the traditional CART and its R implementation `rpart`.

## Overview

Compared to other similar trees, `LDATree` sets itself apart in the
following ways:

- It applies the idea of LDA (Linear Discriminant Analysis) when
  selecting variables, finding splits, and fitting models in terminal
  nodes.

- It addresses certain limitations of the R implementation of LDA
  (`MASS::lda`), such as handling missing values, dealing with more
  features than samples, and constant values within groups.

- Re-implement LDA using the Generalized Singular Value Decomposition
  (GSVD), LDATree offers quick response, particularly with large
  datasets.

- The package also includes several visualization tools to provide
  deeper insights into the data.

## Installation

``` r
install.packages("LDATree")
```

The CRAN version is an outdated one from 08/2023. Please stay tune for
the latest version, which will be released around 10/2024. Meanwhile,
feel free to try the undocumented version bellow.

``` r
library(devtools)
install_github('Moran79/LDATree')
```

## Usage

To build an LDATree:

``` r
library(LDATree)
set.seed(443)
diamonds <- as.data.frame(ggplot2::diamonds)[sample(53940, 2000),]
datX <- diamonds[, -2]
response <- diamonds[, 2] # we try to predict "cut"
fit <- Treee(datX = datX, response = response, verbose = FALSE)
```

To plot the LDATree:

``` r
# View the overall tree.
plot(fit)
```

<img src="man/figures/README-plot1-1.png" width="80%" style="display: block; margin: auto;" />

``` r
# Three types of individual plots
# 1. Scatter plot on first two LD scores
plot(fit, datX = datX, response = response, node = 1)
```

<img src="man/figures/README-plot2-1.png" width="80%" style="display: block; margin: auto;" />

``` r

# 2. Density plot on the first LD score
plot(fit, datX = datX, response = response, node = 3)
```

<img src="man/figures/README-plot2-2.png" width="80%" style="display: block; margin: auto;" />

``` r

# 3. A message
plot(fit, datX = datX, response = response, node = 2)
#> [1] "Every observation in node 2 is predicted to be Fair"
```

To make predictions:

``` r
# Prediction only.
predictions <- predict(fit, datX)
head(predictions)
#> [1] "Ideal" "Ideal" "Ideal" "Ideal" "Ideal" "Ideal"
```

``` r
# A more informative prediction
predictions <- predict(fit, datX, type = "all")
head(predictions)
#>   response node         Fair        Good Very Good     Premium     Ideal
#> 1    Ideal    6 4.362048e-03 0.062196349 0.2601145 0.056664046 0.6166630
#> 2    Ideal    6 1.082022e-04 0.006308281 0.1290079 0.079961227 0.7846144
#> 3    Ideal    6 7.226446e-03 0.077434549 0.2036148 0.023888946 0.6878352
#> 4    Ideal    6 1.695119e-02 0.115233616 0.1551836 0.008302145 0.7043295
#> 5    Ideal    6 4.923729e-05 0.004157352 0.1498265 0.187391975 0.6585749
#> 6    Ideal    6 4.827312e-03 0.061274797 0.1978061 0.027410359 0.7086815
```

## Getting help

If you encounter a clear bug, please file an issue with a minimal
reproducible example on
[GitHub](https://github.com/Moran79/LDATree/issues)
