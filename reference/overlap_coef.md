# Overlapping coefficient between multivariate t-distributions

Compute a Monte Carlo approximation of the overlapping coefficient
between two multivariate t-distributions with arbitrary mean, covariance
and degrees of freedom.

## Usage

``` r
overlap_coef(mean1, mean2, cov1, cov2, df1, df2, nb_sample = 10000)
```

## Arguments

- mean1:

  A vector, the mean parameter of a multi t-distribution

- mean2:

  A vector, the mean parameter of the other multi t-distribution

- cov1:

  A matrix, the covariance parameter of a multi t-distribution

- cov2:

  A matrix, the covariance parameter of the other multi t-distribution

- df1:

  A number, the degrees of freedom of a multi t-distribution

- df2:

  A number, the degrees of freedom of the other multi t-distribution

- nb_sample:

  A number of samples drawn to compute the Monte Carlo estimation

## Value

A number, the Monte Carlo approximation of the overlapping coefficient
between the two multivariate t-distributions.

## Examples

``` r
TRUE
#> [1] TRUE
```
