#' Mulivariate posterior distribution of the means
#'
#' Compute the multivariate posterior distribution of the means
#' between multiple groups, for multiple correlated peptides. The function
#' accounts for multiple imputations through the \code{Draw} identifier in the
#' dataset.
#'
#' @param data A tibble or data frame containing imputed data sets for all
#'    groups. Required columns: \code{Peptide}, \code{Group}, \code{Sample},
#'    \code{Output}. If missing data have been estimated from multiple
#'    imputations, each imputation should be identified in an optional
#'    \code{Draw} column.
#' @param mu_0 A vector, corresponding to the prior mean.
#' @param lambda_0 A number, corresponding to the prior covariance scaling
#'    parameter.
#' @param Sigma_0 A matrix, corresponding to the prior covariance parameter.
#' @param nu_0 A number, corresponding to the prior degrees of freedom.
#'
#' @return A vector providing the empirical multivariate posterior distribution
#'    for the mean of the considered groups.
#' @export
#'
#' @examples
#' TRUE
multi_posterior_mean = function(
    data,
    mu_0 = 0,
    lambda_0 = 1,
    Sigma_0 = NULL,
    nu_0 = 10
){
  t1 = Sys.time()

  ## Loop over the groups
  floop_k = function(k){

    ## Extract the adequate group
    data_k = data %>% dplyr::filter(Group == k)

    ## Collect all the different groups
    list_draw = data_k$Draw %>% unique()

    n_draw = data_k$Draw %>% dplyr::n_distinct()

    list_mat = lapply(list_draw, floop_d, k = k)

    ((1/n_draw) * Reduce('+', list_mat)) %>%
      `colnames<-`(data$Peptide %>% unique()) %>%
      tibble::as_tibble() %>%
      tidyr::pivot_longer(everything(), names_to = 'Peptide', values_to ='Mean') %>%
      dplyr::arrange(Peptide) %>%
      dplyr::mutate(Group = k) %>%
      return()
  }

  ## Loop over the draws
  floop_d = function(d, k){
    t2 = Sys.time()
    paste0('Group n°',k , ' - Draw n° ', d, ' - ', t2 - t1) %>% print()

    ## Extract the adequate draws
    data_k_d = data %>%
      dplyr::filter(Group == k,  Draw == d)

    N_k = data_k_d$Draw %>% dplyr::n_distinct()

    list_draw = data_k_d$Draw %>% unique()

    ## Compute the mean 1/N sum_1^N{y_n}
    mean_yn_k = data_k_d %>%
      dplyr::group_by(Peptide) %>%
      dplyr::summarise(Output = mean(Output)) %>%
      dplyr::pull(Output)

    ##Compute the mean 1/N sum_1^N{(y_n - \bar{y_n}^)2}
    cov_yn = 0
    for(n in list_draw)
    {
      yn_k = data_k_d %>%
        dplyr::filter(Draw == n) %>%
        dplyr::pull(Output)

      centred_y = (yn_k - mean_yn_k)

      cov_yn = cov_yn + tcrossprod(centred_y)
    }

    if(Sigma_0 %>% is.null()){
      Sigma_0 <- diag(1, length(yn_k))
    }

    centred_mean = mean_yn_k - mu_0
    ## Compute the updated posterior hyper-parameters
    mu_N = (lambda_0 * mu_0 + N_k * mean_yn_k) / (lambda_0 + N_k)

    lambda_N = lambda_0 + N_k

    Sigma_N = Sigma_0 + cov_yn +
      ((lambda_0 * N_k) / lambda_N) * tcrossprod(centred_mean)

    nu_N = nu_0 + N_k

    P = length(mu_0) # dimension of the vectors and matrices

    ## Draw from the adequate T-distribution
    # mvtnorm::rmvt(n = n_samples,
    #      sigma = Sigma_N / ((nu_N - P + 1) * lambda_N),
    #      df = nu_N - P + 1,
    #      delta = mu_N) %>%
    #   return()

    tibble::tibble(
      'mu' = mu_N,
      'lambda_N' = lambda_N,
      'Sigma' = Sigma_N,
      'nu' = nu_N
    ) %>% return()

  }

  ## Collect all the different groups
  data$Group %>%
    unique() %>%
    lapply(floop_k) %>%
    dplyr::bind_rows() %>%
    return()
}

#' Posterior distribution of the means
#'
#' Compute the posterior distribution of the means between multiple groups.
#' All peptides are considered independent from one another.
#'
#' @param data A tibble or data frame containing imputed data sets for all
#'   groups. Required columns:  \code{Peptide}, \code{Output}, \code{Group},
#'    \code{Sample}.
#' @param mu_0 A vector, corresponding to the prior mean.
#' @param lambda_0 A number, corresponding to the prior covariance scaling
#'   parameter.
#' @param beta_0 A matrix, corresponding to the prior covariance parameter.
#' @param alpha_0 A number, corresponding to the prior degrees of freedom.
#'
#' @return A vector providing the empirical posterior distribution for the
##   mean of the considered groups.
#' @export
#'
#' @examples
#' TRUE
posterior_mean = function(
    data,
    mu_0 = 0,
    lambda_0 = 1,
    beta_0 = 1,
    alpha_0 = 1
){
  db %>%
    dplyr::group_by(Peptide, Group) %>%
    dplyr::mutate(N_k = dplyr::n_distinct(Sample)) %>%
    dplyr::mutate(SSE = sum( (Output - mean(Output))^2 ) ) %>%
    dplyr::summarise(
      mu = (lambda_0 * mu_0 + N_k * mean(Output)) / (lambda_0 + N_k),
      lambda = lambda_0 + N_k,
      alpha = alpha_0 + (N_k / 2),
      beta = beta_0 + (0.5 * SSE) +
        ((lambda_0 * N_k) / (2 * (lambda_0 + N_k))) * (mean(Output) - mu_0)^2
    ) %>%
    unique() %>%
    return()
}


#' Sample from a t-distribution
#'
#' Sample from a (possibly multivariate) t-distribution. This function can be
#' used to sample both from a prior or posterior, depending on the value of
#' parameters provided.
#'
#' @param posterior A tibble or data frame, detailing for each \code{Peptide}
#'    and each \code{Group}, the value of the t-distribution parameters. The
#'    expected format is typically a return from a posterior_mean() function.
#'    Expected columns in the univariate case: \code{mu}, \code{lambda},
#'    \code{alpha}, \code{beta}. Expected columns in the multivariate case:
#'    \code{mu}, \code{lambda}, \code{Sigma}, \code{nu}.
#' @param nb_sample A number, indicating the number of samples generated for
#'    each couple \code{Peptide}-\code{Group}.
#'
#' @return A tibble containing the \code{Peptide}, \code{Group} and
#'    \code{Sample} columns. The samples of each \code{Peptide}-\code{Group}
#'    couple provide an empirical t-distribution that can be used to compute and
#'    display differences between groups.
#' @export
#'
#' @examples
#' TRUE
sample_distrib = function(posterior, nb_sample = 1000){

  if('Sigma' %in% names(posterior)){
    dist = posterior %>%
      dplyr::group_by(Peptide, Group) %>%
      dplyr::summarise('Sample' = mvtnorm::rmvt(
        n = nb_sample,
        sigma = Sigma / ((nu - P + 1) * lambda),
        df = nu - P + 1,
        delta = mu)
        )
  }

  if('beta' %in% names(posterior)){
    dist = posterior %>%
      dplyr::group_by(Peptide, Group) %>%
      dplyr::summarise('Sample' = mu + sqrt( beta / (lambda * alpha) ) *
                         rt(n = nb_sample, df = 2 * alpha)
      )
  }
    return(dist)
}

