# Vectorised version of multi_posterior_mean()

Alternative vectorised version, highly efficient when nb_peptide \< 30.

## Usage

``` r
vectorised_multi(data, mu_0 = NULL, lambda_0 = 1, Sigma_0 = NULL, nu_0 = 10)
```

## Arguments

- data:

  A tibble or data frame containing imputed data sets for all groups.
  Required columns: `Peptide`, `Group`, `Sample`, `Output`. If missing
  data have been estimated from multiple imputations, each imputation
  should be identified in an optional `Draw` column.

- mu_0:

  A vector, corresponding to the prior mean. If NULL, all groups are
  initialised with the same empirical mean for each peptide.

- lambda_0:

  A number, corresponding to the prior covariance scaling parameter.

- Sigma_0:

  A matrix, corresponding to the prior covariance parameter. If NULL,
  the identity matrix will be used by default.

- nu_0:

  A number, corresponding to the prior degrees of freedom.

## Value

A tibble providing the parameters of the posterior t-distribution for
the mean of the considered groups for each peptide.

## Examples

``` r
TRUE
#> [1] TRUE
```
