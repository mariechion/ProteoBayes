# Identify differences in multivariate posteriors

Compute a multivariate inference criterion to examine whether the
posterior multivariate t-distributions of groups should be considered
different enough to be called 'differential'. Two groups are considered
can be discriminated based on the probability weights for the
element-wise means to be greater in each group. The Credible Intervals
for each marginals are also provided.

## Usage

``` r
multi_identify_diff(
  posterior,
  plot = TRUE,
  overlap_coef = TRUE,
  cumulative = FALSE,
  nb_sample = 1e+05,
  nb_sample_overlap = 10000
)
```

## Arguments

- posterior:

  A tibble, typically coming from a
  [`posterior_mean()`](https://mariechion.github.io/ProteoBayes/reference/posterior_mean.md)
  function, containing the parameters of the multivariate posterior
  t-distributions for the mean of the considered groups and draws for
  each peptide.

- plot:

  A boolean, indicating whether a results plot should be displayed.

- overlap_coef:

  A boolean, indicating whether the overlapping coefficient between
  multivariate t-distributions should be computed for all groups.

- cumulative:

  A boolean, indicating whether the probability distribution should be
  the cumulative instead.

- nb_sample:

  A number of samples to draw from the empirical distributions

- nb_sample_overlap:

  A number of samples to draw for the Monte Carlo approximation of the
  Overlapping Coefficients between all posteriors.

## Value

A list, containing:

- Diff_mean, a tibble containing the posterior mean (and their
  difference) for each peptide and for all one-by-one group comparisons.

- Diff_proba, a tibble containing the probability distribution of the
  number of differential peptides for all one-by-one group comparisons.

- Overlap_coef, a tibble containing the overlapping coefficient between
  the posterior multivariate t-distributions for all one-by-one group
  comparisons.

## Examples

``` r
TRUE
#> [1] TRUE
```
