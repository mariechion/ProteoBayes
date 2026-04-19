# Identify posterior mean differences

Compute a criterion based on Credible Intervals (CI) to determine
whether the posterior t-distributions of groups should be considered
different enough to deserve further examination. Two groups are
considered probably 'distinct' if the Credible Interval of level
`CI_level` of their respective posterior t-distributions do not overlap.

## Usage

``` r
identify_diff(posterior)
```

## Arguments

- posterior:

  A tibble, typically coming from a
  [`posterior_mean()`](https://mariechion.github.io/ProteoBayes/reference/posterior_mean.md)
  function, containing the parameters of the multivariate posterior
  t-distributions for the mean of the considered groups and draws for
  each peptide.

## Value

A tibble, indicating which peptides and groups seem to be different

## Examples

``` r
TRUE
#> [1] TRUE
```
