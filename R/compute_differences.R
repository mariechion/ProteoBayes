#' Identify posterior mean differences
#'
#' Compute a criterion based on Credible Intervals (CI) to determine whether
#' the posterior t-distributions of groups should be considered different enough
#' to deserve further examination. Two groups are considered probably 'distinct'
#' if the Credible Interval of level \code{CI_level} of their respective
#' posterior t-distributions do not overlap.
#'
#' @param posterior A tibble, typically coming from a \code{posterior_mean()}
#'     function, containing the parameters of the multivariate posterior
#'     t-distributions for the mean of the considered groups and draws for each
#'     peptide.
#' @param CI_level A number, defining the order of quantile chosen to assess
#'     differences between groups.
#' @param nb_samples A number (optional), indicating the
#'     number of samples to draw from the posteriors for computing mean and
#'     credible intervals . Only used if \code{posterior} is multivariate,
#'     typically coming from a \code{multi_posterior_mean()} function.
#'
#' @return A tibble, indicating which peptides and groups seem to be different
#' @export
#'
#' @examples
#' TRUE
identify_diff <- function(
    posterior,
    CI_level = 0.05,
    nb_samples = 1000){

  ## We could think of defining a multivariate criterion for differences

  if('Sigma' %in% names(posterior)){
    db_diff <- posterior %>%
      sample_distrib(10) %>%
      dplyr::group_by(.data$Peptide, .data$Group) %>%
      dplyr::reframe('mu' = mean(.data$Sample),
                     'CI_inf' = stats::quantile(.data$Sample, CI_level/2),
                     'CI_sup' = stats::quantile(.data$Sample, 1- CI_level/2))

  } else{
    db_diff <- posterior %>%
      dplyr::mutate('sigma' = sqrt(.data$beta / (.data$lambda * .data$alpha)),
                    'df' =  2 * .data$alpha) %>%
      dplyr::reframe(
        .data$Peptide, .data$Group, .data$mu,
        'CI_inf' = extraDistr::qlst(
          p = CI_level/2,
          df =.data$df,
          mu = .data$mu,
          sigma = .data$sigma,
          lower.tail = T),
        'CI_sup' =  extraDistr::qlst(
          p = CI_level/2,
          df =.data$df,
          mu = .data$mu,
          sigma = .data$sigma,
          lower.tail = F)
      )
  }

    db_diff %>%
      tidyr::expand_grid('Group2' = unique(.data$Group)) %>%
      dplyr::left_join(
        db_diff %>%
          dplyr::rename('Group2' = .data$Group,
                        'mu2' = .data$mu,
                        'CI_inf2' = .data$CI_inf,
                        'CI_sup2' = .data$CI_sup),
        by = c("Peptide", 'Group2')) %>%
      dplyr::filter(.data$Group != .data$Group2) %>%
      dplyr::mutate( 'Distinct' = (.data$CI_sup < .data$CI_inf2) |
                       (.data$CI_inf > .data$CI_sup2)
      ) %>%
      return()

}

