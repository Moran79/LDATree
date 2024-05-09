
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LDATree <a href="http://iamwangsiyu.com/LDATree/"><img src="man/figures/logo.png" align="right" height="139" alt="LDATree website" /></a>

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/LDATree)](https://CRAN.R-project.org/package=LDATree)
[![R-CMD-check](https://github.com/Moran79/LDATree/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Moran79/LDATree/actions/workflows/R-CMD-check.yaml)
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
install.packages("LDATree") # This version is an outdated one from 08/2023.

# As of 06/2024, please use the command below for the current version,
# the official CRAN release will be coming soon!

# library(devtools)
# install_github('Moran79/LDATree')
```

## Usage

To build an LDATree:

``` r
library(LDATree)
set.seed(443)
mpg <- as.data.frame(ggplot2::mpg)
datX <- mpg[, -5] # All predictors without Y
response <- mpg[, 5] # we try to predict "cyl" (number of cylinders)
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
#> [1] "Every observation in this node is predicted to be 4"
```

To make predictions:

``` r
# Prediction only.
predictions <- predict(fit, datX)
head(predictions)
#> [1] "4" "4" "4" "4" "6" "6"
```

``` r
# A more informative prediction
predictions <- predict(fit, datX, type = "all")
head(predictions)
#>   response node 4 5 6 8
#> 1        4   14 1 0 0 0
#> 2        4    6 1 0 0 0
#> 3        4    6 1 0 0 0
#> 4        4    6 1 0 0 0
#> 5        6   18 0 0 1 0
#> 6        6   15 0 0 1 0
```

## Getting help

If you encounter a clear bug, please file an issue with a minimal
reproducible example on
[GitHub](https://github.com/Moran79/LDATree/issues)
