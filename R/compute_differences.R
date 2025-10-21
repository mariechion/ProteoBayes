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
#' @param plot A boolean, indicating whether a results plot should be displayed.
#' @param overlap_coef A boolean, indicating whether the overlapping coefficient
#'     between multivariate t-distributions should be computed for all groups.
#' @param cumulative A boolean, indicating whether the probability distribution
#'     should be the cumulative instead.
#' @param nb_sample A number of samples to draw from the empirical distributions
#' @param nb_sample_overlap A number of samples to draw for the Monte Carlo
#'    approximation of the Overlapping Coefficients between all posteriors.
#'
#' @return A list, containing:
#'    - Diff_mean, a tibble containing the posterior mean (and their difference)
#'    for each peptide and for all one-by-one group comparisons.
#'    - Diff_proba, a tibble containing the probability distribution of the
#'    number of differential peptides for all one-by-one group comparisons.
#'    - Overlap_coef, a tibble containing the overlapping coefficient between
#'    the posterior multivariate t-distributions for all one-by-one group
#'    comparisons.
#'
#' @export
#'
#' @examples
#' TRUE
multi_identify_diff <- function(
    posterior,
    plot = TRUE,
    overlap_coef = TRUE,
    cumulative = FALSE,
    nb_sample = 1e5,
    nb_sample_overlap = 1e4){

  #### Compute the difference of means between groups for each Peptide

  ## Compute the average over all Draws
  averaged_mean = posterior %>%
    dplyr::select(.data$Draw, .data$Group, .data$Peptide, .data$mu) %>%
    dplyr::group_by(.data$Group, .data$Peptide) %>%
    dplyr::summarise(mu = mean(.data$mu), .groups = 'drop')

  ## Compute the difference between means of all Groups
  diff_mean = averaged_mean %>%
    tidyr::expand_grid('Group2' = unique(.data$Group)) %>%
    dplyr::left_join(averaged_mean %>%
                       dplyr::select(.data$Group, .data$mu, .data$Peptide) %>%
                       dplyr::rename(Group2 = .data$Group, mu2 = .data$mu),
                     by = c('Group2', 'Peptide')) %>%
    dplyr::mutate(Diff_mean = .data$mu - .data$mu2) %>%
    dplyr::rename(Mean = .data$mu, Mean2 = .data$mu2) %>%
    dplyr::relocate(.data$Peptide, .before = 1)

  #### Compute a meaningfull multivariate uncertainty measure:
  #### Probability that Group1 has N Peptides greater than in Group2, for all N

  ## Get the dimension
  P = posterior$Peptide %>% unique() %>% length()

  ## Get the number of imputed datasets (Draws)
  nb_draw = posterior$Draw %>% unique() %>% length()

  ## Get the list of Groups
  list_groups = posterior$Group %>% unique() %>% as.character()

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
        dplyr::filter(.data$Group == i) %>%
        dplyr::filter(.data$Draw == j)

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
  diff_proba = tibble::tibble()

  ## Initialise the list of remaining groups to avoid duplicated comparisons
  list_remaining_groups = list_groups

  ## Initialise the tibble of overlap coef for each one-by-one group comparison
  overlap = tibble::tibble()

  ## Loop over all one-by-one group comparisons
  for(i in list_groups){

    ## Remove current group1 from the list of group2
    list_remaining_groups = list_remaining_groups[-1]

    for(j in list_remaining_groups){

      ## Compute the difference between group1 and group2 for each sample
      diff_samples = samples[[i]] - samples[[j]]

      diff_proba = diff_proba %>%
        dplyr::bind_rows(
          diff_samples %>%
            sign() %>%
            (`+`)(1) %>%
            (`/`)(2) %>%
            rowSums() %>%
            tibble::as_tibble() %>%
            dplyr::count(.data$value) %>%
            dplyr::mutate(n = .data$n/nb_sample) %>%
            dplyr::arrange(.data$value) %>%
            dplyr::rename('Nb_peptides' = .data$value, 'Proba' = .data$n) %>%
            tidyr::complete('Nb_peptides'= tidyr::full_seq(0:P, 1),
                            fill= list('Proba' = 0)) %>%
            dplyr::arrange(.data$Nb_peptides) %>%
            dplyr::mutate('Group1' = i) %>%
            dplyr::mutate('Group2' = j)
        )

      if(overlap_coef){
        if(nb_draw > 1){
          cat('Impossible to compute the Overlapping Coefficient when the',
              'number of Draws > 1. There is no mathematical formula to',
              'combine multivariate t-distribution in this case.')
          overlap = NULL
        } else{
          ## Extract the posterior parameters of the corresponding Groups
          db_group1 = posterior %>%
            dplyr::filter(.data$Group == i)

          db_group2 = posterior %>%
            dplyr::filter(.data$Group == j)

          ## Compute posterior degrees of freedom of the multi t-distributions
          df1 = unique(db_group1$nu) - P + 1
          df2 = unique(db_group2$nu) - P + 1

          ## Get the posterior scaling parameters lambda
          lambda1 = unique(db_group1$lambda)
          lambda2 = unique(db_group2$lambda)

          ## Compute posterior covariance matrices of the multi t-distributions
          cov1 = matrix(db_group1$Sigma / (lambda1 * df1), nrow = P, ncol = P)
          cov2 = matrix(db_group2$Sigma / (lambda2 * df2), nrow = P, ncol = P)

          ## Get the posterior mean vector of the multivariate t-distributionS
          mean1 = unique(db_group1$mu)
          mean2 = unique(db_group2$mu)

          tib_overlap = tibble::tibble(
            'Group1' = i,
            'Group2' = j,
            'Overlap_coef' = overlap_coef(
              mean1 = mean1,
              mean2 = mean2,
              cov1 = cov1,
              cov2 = cov2,
              df1 = df1,
              df2 = df2,
              nb_sample = nb_sample_overlap))

          overlap = overlap %>% dplyr::bind_rows(tib_overlap)
        }
      }
    }
  }

  ## Return the cumulative probability distribution if requested
  if(cumulative){
  diff_proba = diff_proba %>%
    dplyr::group_by(.data$Group1, .data$Group2) %>%
    dplyr::mutate(Proba = cumsum(.data$Proba)) %>%
    dplyr::rename('Cumul_proba' = .data$Proba)
  }

  list_diff = list('Diff_mean' = diff_mean, 'Diff_proba' = diff_proba)

  if(overlap_coef){
    list_diff[['Overlap_coef']] = overlap
  }

  ## Display a plot if requested
  if(plot){
    plot_multi_diff(list_diff, cumulative = cumulative)
  }

  list_diff %>%
    return()
}


#' Multidimensional Credible Interval
#'
#' Compute whether a given vector belongs to the multivariate t-distribution
#' Credible Interval for given level, mean, covariance and degrees of freedom.
#'
#' @param data A vector, of compatible dimension with \code{mean} and \code{cov}
#' @param mean A vector, the mean parameter of the multivariate t-distribution
#' @param cov A matrix, the covariance parameter of the multivariate
#'    t-distribution
#' @param df A number, the degrees of freedom of the multivariate t-distribution
#' @param level A number, between 0 and 1, corresponding to the level of the
#'    Credible Interval. Default is 0.95.
#'
#' @returns A boolean, indicating whether the \code{data} vector belongs to the
#'      computed Credible Interval.
#' @export
#'
#' @examples
#' TRUE
multi_CI <- function(
    data,
    mean,
    cov,
    df,
    level = 0.95){

  dim = length(data)

  # Mean term
  m = data - mean

  ## Compute whether 'data' belongs to the Credible Interval
  in_CI = (t(m) %*% solve(cov) %*% m) < dim * stats::qf(level, dim, df)

  return(in_CI)
}

#' Overlapping coefficient between multivariate t-distributions
#'
#' Compute a Monte Carlo approximation of the overlapping coefficient between
#' two multivariate t-distributions with arbitrary mean, covariance and
#' degrees of freedom.
#'
#' @param mean1 A vector, the mean parameter of a multi t-distribution
#' @param mean2 A vector, the mean parameter of the other multi t-distribution
#' @param cov1 A matrix, the covariance parameter of a multi t-distribution
#' @param cov2 A matrix, the covariance parameter of the other multi
#'    t-distribution
#' @param df1 A number, the degrees of freedom of a multi t-distribution
#' @param df2 A number, the degrees of freedom of the other multi t-distribution
#' @param nb_sample A number of samples drawn to compute the Monte Carlo estimation
#'
#' @returns A number, the Monte Carlo approximation of the overlapping
#'    coefficient between the two multivariate t-distributions.
#' @export
#'
#' @examples
#' TRUE
overlap_coef <- function(
    mean1,
    mean2,
    cov1,
    cov2,
    df1,
    df2,
    nb_sample = 1e4
    ) {
    d <- length(mean1)

    # Use an approximate proposal distribution for the importance sampling
    mu_r <- (mean1 + mean2) / 2
    Sigma_r <- (cov1 + cov2) / 2
    nu_r <- min(df1, df2)

    # Generate samples from this approximate proposal distribution
    x_r <- LaplacesDemon::rmvt(n = nb_sample, mu = mu_r, S = Sigma_r, df = nu_r)


    # For numerical stability of small probabilities, log-transform is used
    # Compute the target log-densities (for numerical stability)
    logp <- LaplacesDemon::dmvt(x_r, mu = mean1, S = cov1, df = df1, log = TRUE)
    logq <- LaplacesDemon::dmvt(x_r, mu = mean2, S = cov2, df = df2, log = TRUE)

    # Compute the log-density used to weight the importance sampling
    logr <- LaplacesDemon::dmvt(x_r,
                                mu = mu_r,
                                S = Sigma_r,
                                df = nu_r,
                                log = TRUE)

    # Compute the importance weights (log (min(f,g) / sampling weight) )
    log_weights <- pmin(logp, logq) - logr

    # Compute exponential of the result back to the original scale
    weights <- exp(log_weights)

    # Final Monte Carlo estimate of the overlapping coefficient
    mean(weights) %>%
      return()
}

