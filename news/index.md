# Changelog

## ProteoBayes (development version)

### Major

- Add multi_overlap_coefficient() to compute a Monte Carlo approximation
  of the Overlapping Coefficient in high dimension
- Add a new multi_identify_diff() function for multivariate differential
  analysis
- Update the identify_diff() function to include overlap coefficient
- Add the overlap_coefficient() function to calculate overlap
  coefficient between two sets
- Implement new inference and plotting functions for the multivariate
  results

### Minor

- Fix an issue with the ‘nb_sample’ argument in identify_diff()
- Fix an implementation issue with the mixture of T-distribution in the
  case of multiple Draws

## ProteoBayes 1.0.0

CRAN release: 2023-07-19

- Initial release
