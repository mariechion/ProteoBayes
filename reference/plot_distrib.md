# Plot the posterior distribution of the difference of means

Display the posterior distribution of the difference of means between
two groups for a specific peptide. If only one group is provide, the
function display the posterior distribution of the mean for this
specific group instead. The function provides additional tools to
represent information to help inference regarding the difference between
groups (reference at 0 on the x-axis, probability of `group1` \>
`group2` and conversely).

## Usage

``` r
plot_distrib(
  sample_distrib,
  group1 = NULL,
  group2 = NULL,
  peptide = NULL,
  prob_CI = 0.95,
  show_prob = TRUE,
  mean_bar = TRUE,
  index_group1 = NULL,
  index_group2 = NULL
)
```

## Arguments

- sample_distrib:

  A data frame, typically coming from the
  [`sample_distrib()`](https://mariechion.github.io/ProteoBayes/reference/sample_distrib.md)
  function, containing the following columns: `Peptide`, `Group` and
  `Sample`. This argument should contain the empirical posterior
  distributions to be displayed.

- group1:

  A character string, corresponding to the name of the group for which
  we plot the posterior distribution of the mean. If NULL (default), the
  first group appearing in `sample_distrib` is displayed. If `group2` is
  provided, the posterior difference of the groups is displayed instead.

- group2:

  A character string, corresponding to the name of the group we want to
  compare to `group1`. If NULL (default), only the posterior
  distribution of the mean for `group1` is displayed.

- peptide:

  A character string, corresponding to the name of the peptide for which
  we plot the posterior distribution of the mean. If NULL (default),
  only the first appearing in `sample_distrib` is displayed.

- prob_CI:

  A number, between 0 and 1, corresponding the level of the Credible
  Interval (CI), represented as side regions (in red) of the posterior
  distribution. The default value (0.95) display the 95% CI, meaning
  that the central region (in blue) contains 95% of the probability
  distribution of the mean.

- show_prob:

  A boolean, indicating whether we display the label of probability
  comparisons between the two groups.

- mean_bar:

  A boolean, indicating whether we display the vertical bar
  corresponding to 0 on the x-axis (when comparing two groups), of the
  mean value of the distribution (when displaying a unique group).

- index_group1:

  A character string, used as the index of `group1` in the legends. If
  NULL (default), `group1` is used.

- index_group2:

  A character string, used as the index of `group2` in the legends. If
  NULL (default), `group2` is used.

## Value

Plot of the required posterior distribution.

## Examples

``` r
TRUE
#> [1] TRUE
```
