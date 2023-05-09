#' Generate a synthetic dataset tailored for ProteoBayes
#'
#' Simulate a complete training dataset, which may be representative of various
#' applications. Several flexible arguments allow adjustment of the number of
#' peptides, of groups, and samples in each experiment. The values of several
#' parameters controlling the data generation process can be modified.
#'
#' @param nb_peptide An integer, indicating the number of peptides in the data.
#' @param nb_group An integer, indicating the number of groups/conditions.
#' @param nb_sample  An integer, indicating the number of samples in the data
#'   for each peptide (i.e the repetitions of the same experiment).
#' @param multi_imp A boolean, indicating whether multiple imputations have
#'   been applied to obtain the dataset.
#' @param nb_draw A number, indicating the number of imputation procedures
#'   applied to obtain this dataset.
#' @param range_peptide A 2-sized vector, indicating the range of values from
#'   which to pick a mean value for each peptide.
#' @param diff_group A number, indicating the mean difference between
#'   consecutive groups
#' @param var_sample A number, indicating the noise variance for each new
#'   sample of a peptide.
#' @param var_draw A number, indicating the noise variance for each
#'   imputation draw.
#'
#' @return  A full dataset of synthetic data.
#' @export
#'
#' @examples
#'
#' ## Generate a dataset with 5 peptides in each of the 2 groups, observed for
#' ##  3 different samples
#' data = simu_db(nb_peptide = 5, nb_group = 2, nb_sample = 3)
#'
#' ## Generate a dataset with 3 peptides in each of the 3 groups, observed for
#' ## 4 different samples, for which 5 imputation draws are available.
#' data = simu_db(nb_peptide = 3, nb_group = 3, nb_sample = 4; nb_draws = 5)
#'
simu_db = function(
    nb_peptide = 5,
    nb_group = 2,
    nb_sample = 5,
    multi_imp = FALSE,
    nb_draw = 5,
    range_peptide = c(0, 50),
    diff_group = 3,
    var_sample = 2,
    var_draw = 1
){
  db <- tibble::tibble(
    'Peptide' = rep(paste0('Peptide_', 1:nb_peptide), each= nb_group*nb_sample),
    'Group' = rep(rep(1:nb_group, each = nb_sample),nb_peptide),
    'Sample' = rep(rep( 1:nb_sample, nb_group*nb_peptide))
  ) %>%
    dplyr::group_by(Peptide) %>%
    dplyr::mutate(.data$Output = runif(1,
                                       range_peptide[1],
                                       range_peptide[2])) %>%
    dplyr::group_by(.data$Group) %>%
    dplyr::mutate(.data$Output = .data$Output + diff_group * .data$Group) %>%
    dplyr::group_by(.data$Sample) %>%
    dplyr::mutate(.data$Output = .data$Output + rnorm(1, 0, var_sample)) %>%
    dplyr::ungroup()

  if(multi_imp){
    db = db %>%
      tidyr::uncount(nb_draw, .id = 'Draw') %>%
      dplyr::group_by(.data$Draw) %>%
      dplyr::mutate(.data$Output = .data$Output + rnorm(1, 0, var_draw)) %>%
      dplyr::ungroup()
  }
    db %>% return()
}
