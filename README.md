# bacon

![](https://github.com/evanjflack/bacon/workflows/R-CMD-check/badge.svg)
[![Travis build status](https://travis-ci.org/evanjflack/bacon.svg?branch=master)](https://travis-ci.org/evanjflack/bacon)
[![Coverage status](https://codecov.io/gh/evanjflack/bacon/branch/master/graph/badge.svg)](https://codecov.io/github/evanjflack/bacon?branch=master)
[![Example Jupyter Notebook](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/evanjflack/bacon/master?filepath=index.ipynb)

## Overview

`bacon` is a package with tools to perform the Goodman-Bacon decomposition for differences-in-differences with variation in treatment timing.

## Functions
* `bacon()`: calculates all 2x2 differences-in-differences estimates and weights for the Bacon-Goodman decomposition.

## Data
* `math_refom`: A data set containing state/year level data on an educational reform and future income.
* `castle`: A data set containing state/year level data on a crime reform.

## Installation
``` r
library(devtools)
install_github("evanjflack/bacon")
```

## Usage
```r
df_bacon <- bacon(incearn_ln ~ reform_math,
                  data = bacon::math_reform,
                  id_var = "state",
                  time_var = "class")

library(ggplot2)

ggplot(df_bacon) +
  aes(x = weight, y = estimate, shape = factor(type)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  labs(x = "Weight", y = "Estimate", shape = "Type")
```

## References
Goodman-Bacon, Andrew. 2018. "Difference-in-Differences with Variation in Treatment Timing." National Bureau
of Economic Research Working Paper Series No. 25018. doi: 10.3386/w25018.

[Paper Link](https://cdn.vanderbilt.edu/vu-my/wp-content/uploads/sites/2318/2019/07/29170757/ddtiming_7_29_2019.pdf)
