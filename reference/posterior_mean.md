# Posterior distribution of the means

Compute the posterior distribution of the means between multiple groups.
All peptides are considered independent from one another.

## Usage

``` r
posterior_mean(data, mu_0 = NULL, lambda_0 = 1, beta_0 = 1, alpha_0 = 1)
```

## Arguments

- data:

  A tibble or data frame containing imputed data sets for all groups.
  Required columns: `Peptide`, `Output`, `Group`, `Sample`.

- mu_0:

  A vector, corresponding to the prior mean.

- lambda_0:

  A number, corresponding to the prior covariance scaling parameter.

- beta_0:

  A matrix, corresponding to the prior covariance parameter.

- alpha_0:

  A number, corresponding to the prior degrees of freedom.

## Value

A tibble providing the empirical posterior distribution for the

## Examples

``` r
TRUE
#> [1] TRUE
```
