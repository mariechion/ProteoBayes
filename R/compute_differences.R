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
#' @param nb_sample A number (optional), indicating the
#'     number of samples to draw from the posteriors for computing mean and
#'     credible intervals . Only used if \code{posterior} is multivariate,
#'     typically coming from a \code{multi_posterior_mean()} function.
#' @param separate_groups A boolean, indicating whether the distributions of
#'     groups should be presented separately or directly as a difference.
#'
#' @return A tibble, indicating which peptides and groups seem to be different
#' @export
#'
#' @examples
#' TRUE
identify_diff <- function(
    posterior,
    CI_level = 0.05,
    nb_sample = 1000,
    separate_groups = TRUE){

  if(separate_groups){
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

    db_diff %>%
      tidyr::expand_grid('Group2' = unique(.data$Group)) %>%
      dplyr::left_join(
        db_diff %>%
          dplyr::rename('Group2' = .data$Group,
                        'mu2' = .data$mu,
                        'CI_inf2' = .data$CI_inf,
                        'CI_sup2' = .data$CI_sup),
        by = c("Peptide", 'Group2')) %>%
      dplyr::filter(.data$Group != .data$Group2)
  } else{
    db_diff <- posterior %>%
      dplyr::mutate('sigma' = sqrt(.data$beta / (.data$lambda * .data$alpha)),
                    'df' =  2 * .data$alpha)

    db_diff %>%
      tidyr::expand_grid('Group2' = unique(.data$Group)) %>%
      dplyr::left_join(
        db_diff %>%
          dplyr::rename('Group2' = .data$Group,
                        'mu2' = .data$mu,
                        'CI_inf2' = .data$CI_inf,
                        'CI_sup2' = .data$CI_sup),
        by = c("Peptide", 'Group2')) %>%
      dplyr::filter(.data$Group != .data$Group2)
  }

  return(db_diff)
}

#' Identify differences in multivariate posteriors
#'
#' Compute a multivariate inference criterion to examine whether the posterior
#' multivariate t-distributions of groups should be considered different enough
#' to be called 'differential'. Two groups are considered can be discriminated
#' based on the probability weights for the element-wise means to be greater in
#' each group. The Credible Intervals for each marginals are also provided.
#'
#' @param posterior A tibble, typically coming from a \code{posterior_mean()}
#'     function, containing the parameters of the multivariate posterior
#'     t-distributions for the mean of the considered groups and draws for each
#'     peptide.
#' @param CI_level A number, defining the order of quantile chosen to assess
#'     differences between groups.
#' @param nb_sample A number (optional), indicating the
#'     number of samples to draw from the posteriors for computing mean and
#'     credible intervals . Only used if \code{posterior} is multivariate,
#'     typically coming from a \code{multi_posterior_mean()} function.
#' @param posterior A tibble, typically coming from a \code{posterior_mean()}
#'     function, containing the parameters of the multivariate posterior
#'     t-distributions for the mean of the considered groups and draws for each
#'     peptide.
#' @param CI_level A number, defining the order of quantile chosen to assess
#'     differences between groups.
#' @param nb_sample A number (optional), indicating the
#'     number of samples to draw from the posteriors for computing mean and
#'     credible intervals . Only used if \code{posterior} is multivariate,
#'     typically coming from a \code{multi_posterior_mean()} function.
#' @param plot A boolean, indicating whether a results plot should be displayed.
#' @param cumulative A boolean, indicating whether the probability distribution
#'     should be the cumulative instead.
#' @param nb_sample A number of samples to draw from the empirical distributions
#' @param nb_sample_per_dim
#'
#' @return A tibble, indicating the probability distribution for the number of
#'     differentiable peptides for all one-by-one group comparisons
#' @export
#'
#' @examples
#' TRUE
multi_identify_diff <- function(
    posterior,
    plot = TRUE,
    cumulative = FALSE,
    CI_level = 0.05,
    nb_sample = 10000,
    nb_sample_per_dim = NULL){

  ## Get the dimension
  P = posterior$Peptide %>% unique() %>% length()

  ## Get the number of imputed datasets (Draws)
  nb_draw = posterior$Draw %>% unique() %>% length()

  ## Get the list of Groups
  list_groups = posterior$Group %>% unique()

  ## Throw an error if nu < P - 1
  if(min(posterior$nu) < P - 1){
    stop("The 'nu' parameter is too small compared to the number of Peptides",
         " (nu < P - 1). Consider changing prior value for 'nu' or decreasing",
         " the number of Peptides.")
  }

  ## Initialise the list of samples for all groups
  samples = list()

  ## Loop over all groups to generate samples
  for(i in list_groups){

    ## Initialise the object that will average samples of all imputation draws
    samples_draw = 0

    ## Loop over all imputed datasets to generate samples
    for(j in unique(posterior$Draw)){

      ## Extract the posterior parameters of one group
      db_group = posterior %>%
        dplyr::filter(Group == i) %>%
        dplyr::filter(Draw == j)

      ## Compute the posterior degrees of freedom of the multi t-distribution
      df = unique(db_group$nu) - P + 1

      ## Get the posterior scaling parameter lambda
      lambda = unique(db_group$lambda)

      ## Compute the posterior covariance matrix of the multi t-distribution
      Sigma = matrix(db_group$Sigma / (lambda * df), nrow = P, ncol = P)

      ## Get the posterior mean vector of the multivariate t-distribution
      mu = unique(db_group$mu)

      ## Average samples across all imputed datasets
      samples_draw = samples_draw +
        (1/nb_draw) * mvtnorm::rmvt(n=nb_sample, sigma=Sigma, df=df, delta=mu)
    }
    samples[[i]] = samples_draw
  }

  ## Initialise the tibble of results for each one-by-one group comparison
  res = tibble::tibble()

  ## Initialise the list of remaining groups to avoid duplicated comparisons
  list_remaining_groups = list_groups

  ## Loop over all one-by-one group comparisons
  for(i in list_groups){

    ## Remove current group1 from the list of group2
    list_remaining_groups = list_remaining_groups[-1]

    for(j in list_remaining_groups){

      ## Compute the difference between group1 and group2 for each sample
      diff_samples = samples[[i]] - samples[[j]]

      res = res %>%
        dplyr::bind_rows(
          diff_samples %>%
            sign() %>%
            (`+`)(1) %>%
            (`/`)(2) %>%
            rowSums() %>%
            tibble::as_tibble() %>%
            dplyr::count(value) %>%
            dplyr::mutate(n = n/nb_sample) %>%
            dplyr::arrange(value) %>%
            dplyr::rename('Nb_peptides' = value, 'Proba' = n) %>%
            tidyr::complete('Nb_peptides'= tidyr::full_seq(0:P, 1),
                            fill= list('Proba' = 0)) %>%
            dplyr::arrange(Nb_peptides) %>%
            dplyr::mutate('Group1' = i) %>%
            dplyr::mutate('Group2' = j)
        )
    }
  }

  ## Return the cumulative probability distribution if requested
  if(cumulative){
  res = res %>%
    dplyr::group_by(Group1, Group2) %>%
    dplyr::mutate(Proba = cumsum(Proba)) %>%
    dplyr::rename('Cumul_proba' = Proba)
  }

  ## Display a plot if requested
  if(plot){
    plot_multi_diff(res, cumulative = cumulative)
  }

  res %>%
    return()
}
