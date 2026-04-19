# Overlapping coefficient between univariate t-distributions

Compute a (high speed) quadrature approximation of the overlapping
coefficient between two univariate t-distributions with arbitrary mean,
covariance and degrees of freedom.

## Usage

``` r
overlap_coef(mean1, mean2, var1, var2, df1, df2)
```

## Arguments

- mean1:

  A vector, the mean parameter of a t-distribution

- mean2:

  A vector, the mean parameter of the other t-distribution

- var1:

  A matrix, the variance parameter of a t-distribution

- var2:

  A matrix, the variance parameter of the other t-distribution

- df1:

  A number, the degrees of freedom of a t-distribution

- df2:

  A number, the degrees of freedom of the other t-distribution

- nb_sample:

  A number of samples drawn to compute the Monte Carlo estimation

## Value

A number, the Monte Carlo approximation of the overlapping coefficient
between the two univariate t-distributions.

## Examples

``` r
TRUE
#> [1] TRUE
```
