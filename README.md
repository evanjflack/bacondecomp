
![](https://github.com/evanjflack/bacon/workflows/R-CMD-check/badge.svg) [![Travis build status](https://travis-ci.org/evanjflack/bacon.svg?branch=master)](https://travis-ci.org/evanjflack/bacon) [![Coverage status](https://codecov.io/gh/evanjflack/bacon/branch/master/graph/badge.svg)](https://codecov.io/github/evanjflack/bacon?branch=master) [![Example Jupyter Notebook](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/evanjflack/bacon/master?filepath=index.ipynb)

<!-- README.md is generated from README.Rmd. Please edit that file -->
bacon
=====

`bacon` is a package with tools to perform the Goodman-Bacon decomposition for differences-in-differences with variation in treatment timing. The decomposition can be done with and without time-varying covariates.

Installation
------------

You can install the released version of bacon from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("bacon")
```

Functions
---------

-   `bacon()`: calculates all 2x2 differences-in-differences estimates and weights for the Bacon-Goodman decomposition.

Data
----

-   `math_refom`: Aggregated data from Goodman (2019, JOLE)
-   `castle`: Data from Cheng and Hoekstra (2013, JHR)
-   `divorce:` Data from Stevenson and Wolfers (2006, QJE)

Example
-------

This is a basic example which shows you how to use the bacon() function to decompose the two-way fixed effects estimate of the effect of an education reform on future earnings following Goodman (2019, JOLE).

``` r
library(bacon)
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

<img src="man/figures/README-example-1.png" width="100%" />

References
----------

Goodman-Bacon, Andrew. 2018. "Difference-in-Differences with Variation in Treatment Timing." National Bureau of Economic Research Working Paper Series No. 25018. doi: 10.3386/w25018.

[Paper Link](https://cdn.vanderbilt.edu/vu-my/wp-content/uploads/sites/2318/2019/07/29170757/ddtiming_7_29_2019.pdf)
