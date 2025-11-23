# Multivariate posterior distribution of the means

Compute the multivariate posterior distribution of the means between
multiple groups, for multiple correlated peptides. The function accounts
for multiple imputations through the `Draw` identifier in the dataset.

## Usage

``` r
multi_posterior_mean(
  data,
  mu_0 = NULL,
  lambda_0 = 1,
  Sigma_0 = NULL,
  nu_0 = 10,
  vectorised = FALSE
)
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

- vectorised:

  A boolean, indicating whether we should used a vectorised version of
  the function. Default when nb_peptides \< 30. If nb_peptides \> 30,
  there is a high risk that the vectorised version would be slower.

## Value

A tibble providing the parameters of the multivariate posterior
t-distribution for the mean of the considered groups and draws for each
peptide.

## Examples

``` r
TRUE
#> [1] TRUE
```
