# Sample from a t-distribution

Sample from a (possibly multivariate) t-distribution. This function can
be used to sample both from a prior or posterior, depending on the value
of parameters provided.

## Usage

``` r
sample_distrib(posterior, nb_sample = 1000)
```

## Arguments

- posterior:

  A tibble or data frame, detailing for each `Peptide` and each `Group`,
  the value of the t-distribution parameters. The expected format is
  typically a return from a posterior_mean() function. Expected columns
  in the univariate case: `mu`, `lambda`, `alpha`, `beta`. Expected
  columns in the multivariate case: `mu`, `lambda`, `Sigma`, `nu`.

- nb_sample:

  A number, indicating the number of samples generated for each couple
  `Peptide`-`Group`.

## Value

A tibble containing the `Peptide`, `Group` and `Sample` columns. The
samples of each `Peptide`-`Group` couple provide an empirical
t-distribution that can be used to compute and display differences
between groups.

## Examples

``` r
TRUE
#> [1] TRUE
```
