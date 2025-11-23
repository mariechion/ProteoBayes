# Identify posterior mean differences

Compute a criterion based on Credible Intervals (CI) to determine
whether the posterior t-distributions of groups should be considered
different enough to deserve further examination. Two groups are
considered probably 'distinct' if the Credible Interval of level
`CI_level` of their respective posterior t-distributions do not overlap.

## Usage

``` r
identify_diff(
  posterior,
  CI_level = 0.05,
  nb_sample = 1000,
  separate_groups = TRUE
)
```

## Arguments

- posterior:

  A tibble, typically coming from a
  [`posterior_mean()`](https://mariechion.github.io/ProteoBayes/reference/posterior_mean.md)
  function, containing the parameters of the multivariate posterior
  t-distributions for the mean of the considered groups and draws for
  each peptide.

- CI_level:

  A number, defining the order of quantile chosen to assess differences
  between groups.

- nb_sample:

  A number (optional), indicating the number of samples to draw from the
  posteriors for computing mean and credible intervals . Only used if
  `posterior` is multivariate, typically coming from a
  [`multi_posterior_mean()`](https://mariechion.github.io/ProteoBayes/reference/multi_posterior_mean.md)
  function.

- separate_groups:

  A boolean, indicating whether the distributions of groups should be
  presented separately or directly as a difference.

## Value

A tibble, indicating which peptides and groups seem to be different

## Examples

``` r
TRUE
#> [1] TRUE
```
