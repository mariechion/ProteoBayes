# Multidimensional Credible Interval

Compute whether a given vector belongs to the multivariate
t-distribution Credible Interval for given level, mean, covariance and
degrees of freedom.

## Usage

``` r
multi_CI(data, mean, cov, df, level = 0.95)
```

## Arguments

- data:

  A vector, of compatible dimension with `mean` and `cov`

- mean:

  A vector, the mean parameter of the multivariate t-distribution

- cov:

  A matrix, the covariance parameter of the multivariate t-distribution

- df:

  A number, the degrees of freedom of the multivariate t-distribution

- level:

  A number, between 0 and 1, corresponding to the level of the Credible
  Interval. Default is 0.95.

## Value

A boolean, indicating whether the `data` vector belongs to the computed
Credible Interval.

## Examples

``` r
TRUE
#> [1] TRUE
```
