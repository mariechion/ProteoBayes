# Generate a synthetic dataset tailored for ProteoBayes

Simulate a complete training dataset, which may be representative of
various applications. Several flexible arguments allow adjustment of the
number of peptides, of groups, and samples in each experiment. The
values of several parameters controlling the data generation process can
be modified.

## Usage

``` r
simu_db(
  nb_peptide = 5,
  nb_group = 2,
  nb_sample = 5,
  multi_imp = FALSE,
  nb_draw = 5,
  range_peptide = c(0, 50),
  diff_group = 3,
  var_sample = 2,
  var_draw = 1
)
```

## Arguments

- nb_peptide:

  An integer, indicating the number of peptides in the data.

- nb_group:

  An integer, indicating the number of groups/conditions.

- nb_sample:

  An integer, indicating the number of samples in the data for each
  peptide (i.e the repetitions of the same experiment).

- multi_imp:

  A boolean, indicating whether multiple imputations have been applied
  to obtain the dataset.

- nb_draw:

  A number, indicating the number of imputation procedures applied to
  obtain this dataset.

- range_peptide:

  A 2-sized vector, indicating the range of values from which to pick a
  mean value for each peptide.

- diff_group:

  A number, indicating the mean difference between consecutive groups

- var_sample:

  A number, indicating the noise variance for each new sample of a
  peptide.

- var_draw:

  A number, indicating the noise variance for each imputation draw.

## Value

A full dataset of synthetic data.

## Examples

``` r
## Generate a dataset with 5 peptides in each of the 2 groups, observed for
##  3 different samples
data = simu_db(nb_peptide = 5, nb_group = 2, nb_sample = 3)

## Generate a dataset with 3 peptides in each of the 3 groups, observed for
## 4 different samples, for which 5 imputation draws are available.
data = simu_db(nb_peptide = 3, nb_group = 3, nb_sample = 4, nb_draw = 5)
```
