# Plot multivariate comparison between groups

Graphical representation of inference in a multivariate difference
analysis context. The plotted empirical distribution represents, for
each one-to-one group comparison, the probability of the number of
elements for which the mean of a peptide is larger in a given group.
This provides visual evidence on whether two groups are differential or
not, with adequate uncertainty quantification.

## Usage

``` r
plot_multi_diff(multi_diff, plot_mean = TRUE, cumulative = FALSE)
```

## Arguments

- multi_diff:

  A tibble, typically coming from the `multi_identify_diff` function,
  containing probability distribution (or the cumulative distribution)
  of the number of larger Peptides in each one-to-one group comparisons.

- plot_mean:

  A boolean, indicating whether the graph for difference of means
  between all groups for each Peptide should be displayed.

- cumulative:

  A boolean, indicating whether the cumulative distribution or the
  original probability distribution should be displayed.

## Value

A graph (or a matrix of graphs) displaying the multivariate differential
inference summary between groups

## Examples

``` r
TRUE
#> [1] TRUE
```
