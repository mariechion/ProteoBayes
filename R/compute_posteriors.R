#' Multivariate posterior distribution of the means
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
#' @param mu_0 A vector, corresponding to the prior mean. If NULL, all groups
#'    are initialised with the same empirical mean for each peptide.
#' @param lambda_0 A number, corresponding to the prior covariance scaling
#'    parameter.
#' @param Sigma_0 A matrix, corresponding to the prior covariance parameter.
#'    If NULL, the identity matrix will be used by default.
#' @param nu_0 A number, corresponding to the prior degrees of freedom.
#' @param vectorised A boolean, indicating whether we should used a vectorised
#'    version of the function. Default when nb_peptides < 30.
#'    If nb_peptides > 30, there is a high risk that the vectorised version
#'    would be slower.
#'
#' @return A tibble providing the parameters of the multivariate posterior
#'     t-distribution for the mean of the considered groups and draws for
#'     each peptide.
#' @export
#'
#' @examples
#' TRUE
multi_posterior_mean = function(
    data,
    mu_0 = NULL,
    lambda_0 = 1,
    Sigma_0 = NULL,
    nu_0 = 10,
    vectorised = FALSE
){

  ## Remove missing data if present
  data = data %>%
    drop_na()

  ## Initialise prior mean \mu_0 with empirical mean across all groups
  if(mu_0 %>% is.null()){
    data <- data %>%
      dplyr::group_by(.data$Peptide) %>%
      dplyr::mutate('mu_0' = mean(.data$Output)) %>%
      dplyr::ungroup()
  } else {
    data <- data %>%
      dplyr::mutate('mu_0' = mu_0)
  }


  ## If no 'Draw' column, initialise a dummy one
  if(!('Draw' %in% names(data))){
    data <- data %>%
      dplyr::mutate('Draw' = 1)
  }

  ## Extract the dimension
  dim = data$Peptide %>% unique() %>% length()

  ## If NULL, define Sigma_0 as the identity matrix
  if(Sigma_0 %>% is.null()){
    Sigma_0 <- diag(1, dim)
  }

  ## Use vectorised version (default if nb_peptides < 30)
  if(vectorised){
    return(
      vectorised_multi(
        data = data,
        mu_0 = mu_0,
        lambda_0,
        Sigma_0 = Sigma_0,
        nu_0 = nu_0) %>%
      dplyr::arrange(.data$Group)
    )
  }

  ## Loop over the groups
  floop_k = function(k){

    ## Extract the adequate group
    data_k = data %>% dplyr::filter(.data$Group == k)

    ## Collect all the different groups
    list_draw = data_k$Draw %>% unique()

    n_draw = data_k$Draw %>% dplyr::n_distinct()

    lapply(list_draw, floop_d, k = k) %>%
      dplyr::bind_rows() %>%
      tidyr::pivot_longer(!c('Draw', 'Peptide', 'mu', 'lambda','nu'),
                          names_to = 'Peptide2', values_to = 'Sigma') %>%
      dplyr::relocate(.data$Peptide2, .after = .data$Peptide) %>%
      dplyr::mutate("Group" = k, .after = .data$Draw) %>%
      return()
  }

  ## Loop over the draws
  floop_d = function(d, k){
    ## Extract the adequate draws
    data_k_d = data %>%
      dplyr::filter(.data$Group == k,  .data$Draw == d) %>%
      dplyr::arrange(.data$Peptide)

    ## Extract the list of Peptide
    list_Peptide = data_k_d %>%
      dplyr::pull(.data$Peptide) %>%
      unique()

    N_k = data_k_d$Sample %>% dplyr::n_distinct()
    list_sample = data_k_d$Sample %>% unique()

    ## Compute the mean 1/N sum_1^N{y_n}
    mean_yn_k = data_k_d %>%
      dplyr::group_by(.data$Peptide) %>%
      dplyr::reframe("C_Output" = mean(.data$Output),
                     "mu_0" = unique(.data$mu_0))


    ## Compute sum_1^N{(y_n - \bar{y_n}^)2}
    cov_yn = 0
    for(n in list_sample)
    {
      yn_k = data_k_d %>%
        dplyr::filter(.data$Sample == n) %>%
        dplyr::pull(.data$Output)

      centred_y = (yn_k - mean_yn_k$C_Output)

      cov_yn = cov_yn + tcrossprod(centred_y)
    }

    centred_mean = mean_yn_k %>%
      dplyr::reframe('c_mean' = .data$C_Output - .data$mu_0) %>%
      dplyr::pull(.data$c_mean)

    ## Compute the updated posterior hyper-parameters
    mu_N <- (lambda_0*mean_yn_k$mu_0 + N_k*mean_yn_k$C_Output)/(lambda_0 + N_k)

    lambda_N <- lambda_0 + N_k

    nu_N <- nu_0 + N_k

    Sigma_N <- Sigma_0 + cov_yn +
      ((lambda_0 * N_k) / lambda_N) * tcrossprod(centred_mean)
    colnames(Sigma_N) <- list_Peptide

    tibble::tibble(
      'Draw' = d,
      'Peptide' = list_Peptide,
      'mu' = mu_N,
      'lambda' = lambda_N,
      'nu' = nu_N,
      Sigma_N %>% tibble::as_tibble()
    ) %>% return()

  }
  ## Collect all the different groups
  data$Group %>%
    unique() %>%
    lapply(floop_k) %>%
    dplyr::bind_rows() %>%
    return()
}

#' Vectorised version of multi_posterior_mean()
#'
#' Alternative vectorised version, highly efficient when nb_peptide < 30.
#'
#' @param data A tibble or data frame containing imputed data sets for all
#'    groups. Required columns: \code{Peptide}, \code{Group}, \code{Sample},
#'    \code{Output}. If missing data have been estimated from multiple
#'    imputations, each imputation should be identified in an optional
#'    \code{Draw} column.
#' @param mu_0 A vector, corresponding to the prior mean. If NULL, all groups
#'    are initialised with the same empirical mean for each peptide.
#' @param lambda_0 A number, corresponding to the prior covariance scaling
#'    parameter.
#' @param Sigma_0 A matrix, corresponding to the prior covariance parameter.
#'    If NULL, the identity matrix will be used by default.
#' @param nu_0 A number, corresponding to the prior degrees of freedom.
#'
#' @return A tibble providing the parameters of the posterior t-distribution for
#'    the mean of the considered groups for each peptide.
#'
#' @examples
#' TRUE
vectorised_multi = function(
    data,
    mu_0 = NULL,
    lambda_0 = 1,
    Sigma_0 = NULL,
    nu_0 = 10
){

  ## Remove missing data if present
  data = data %>%
    drop_na()

  ## Initialise prior mean \mu_0 with empirical mean across all groups
  if(mu_0 %>% is.null()){
    data <- data %>%
      dplyr::group_by(.data$Peptide) %>%
      dplyr::mutate('mu_0' = mean(.data$Output)) %>%
      dplyr::ungroup()
  } else {
    data <- data %>%
      dplyr::mutate('mu_0' = mu_0)
  }

  ## If no 'Draw' column, initialise a dummy one
  if(!('Draw' %in% names(data))){
    data <- data %>%
      dplyr::mutate('Draw' = 1)
  }

  ## Create a vectorised version of the Sigma_0 matrix
  db_Sigma_0 <- Sigma_0 %>%
    `colnames<-`(unique(data$Peptide)) %>%
    `rownames<-`(unique(data$Peptide)) %>%
    tibble::as_tibble(rownames = 'Peptide') %>%
    tidyr::pivot_longer(- 'Peptide',
                        names_to = 'Peptide2',
                        values_to = 'Sigma_0')

  ## Compute the centred vectors and centred means across all Samples
  mat <- data %>%
    dplyr::group_by(.data$Draw, .data$Group, .data$Peptide) %>%
    dplyr::mutate('N_k' = dplyr::n_distinct(.data$Sample)) %>%
    dplyr::mutate('Mean_Output' = mean(.data$Output)) %>%
    dplyr::mutate('C_Output' = .data$Output - .data$Mean_Output) %>%
    dplyr::mutate('C_Mean' = .data$Mean_Output - .data$mu_0)

  ## Compute the univariate posterior parameters
  post <- mat %>%
    dplyr::reframe(
      'mu' = (lambda_0 * .data$mu_0 + .data$N_k * .data$Mean_Output) /
        (lambda_0 +.data$N_k),
      'lambda' = lambda_0 + .data$N_k,
      'nu' = nu_0 + .data$N_k,
      .data$N_k
    ) %>%
    dplyr::distinct()

  ## Compute the posterior matrix parameter Sigma_N and bind with the others
  mat %>%
    tidyr::expand_grid('Peptide2' = unique(.data$Peptide)) %>%
    dplyr::left_join(
      mat %>%
        dplyr::select(- c('Output', 'mu_0', 'N_k', 'Mean_Output')) %>%
        dplyr::rename('Peptide2' = .data$Peptide,
                      'C_Output2' = .data$C_Output,
                      'C_Mean2' = .data$C_Mean),
      by = c("Group", "Sample", "Draw", "Peptide2")) %>%
    dplyr::mutate('C_i' = .data$C_Output * .data$C_Output2,
                  'C_0' = .data$C_Mean * .data$C_Mean2) %>%
    dplyr::group_by(.data$Draw, .data$Group, .data$Peptide, .data$Peptide2) %>%
    dplyr::reframe('C_i' = sum(.data$C_i),
                   'C_0' = mean(.data$C_0)
                  ) %>%
    dplyr::left_join(db_Sigma_0, by = c("Peptide", "Peptide2")) %>%
    dplyr::left_join(post, by = c("Draw", "Group", "Peptide")) %>%
    dplyr::reframe(.data$Draw, .data$Group, .data$Peptide, .data$Peptide2,
                   .data$mu, .data$lambda, .data$nu,
                   'Sigma' = .data$Sigma_0 + .data$C_i +
                       ((.data$N_k * lambda_0) / .data$lambda) * .data$C_0,
                     ) %>%
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
#' @return A tibble providing the empirical posterior distribution for the
##   mean of the considered groups.
#' @export
#'
#' @examples
#' TRUE
posterior_mean = function(
    data,
    mu_0 = NULL,
    lambda_0 = 1,
    beta_0 = 1,
    alpha_0 = 1
){

  ## Remove missing data if present
  data = data %>%
    drop_na()

  ## Initialise prior mean \mu_0 with empirical mean across all groups
  if(mu_0 %>% is.null()){
    data <- data %>%
      dplyr::group_by(.data$Peptide) %>%
      dplyr::mutate('mu_0' = mean(.data$Output)) %>%
      dplyr::ungroup()
  } else {
    data <- data %>%
      dplyr::mutate('mu_0' = mu_0)
  }

 data %>%
    dplyr::group_by(.data$Peptide, .data$Group) %>%
    dplyr::reframe('N_k' = dplyr::n_distinct(.data$Sample),
                   'mean_output' = mean(.data$Output),
                   'C_i' = sum( (.data$Output - .data$mean_output)^2 ),
                   'mu_0' = unique(.data$mu_0)) %>%
    dplyr::reframe(
      .data$Peptide,
      .data$Group,
      mu = (lambda_0 * .data$mu_0 + .data$N_k*.data$mean_output) /
        (lambda_0 + .data$N_k),
      lambda = lambda_0 + .data$N_k,
      alpha = alpha_0 + (.data$N_k / 2),
      beta = beta_0 + (0.5 * .data$C_i) +
        ((lambda_0 * .data$N_k) / (2 * (lambda_0 + .data$N_k))) *
        (.data$mean_output - .data$mu_0)^2
    ) %>%
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

  if('nu' %in% names(posterior)){

    ## Extract the dimension
    P = posterior$Peptide %>% unique() %>% length()

    ## Throw an error if nu < P - 1
    if(min(posterior$nu) < P - 1){
      stop("The 'nu' parameter is too small compared to the number of Peptides",
      " (nu < P - 1). Consider increasing the prior value for 'nu' or ",
      " decreasing the number of Peptides.")
    }

    dist = posterior %>%
      dplyr::mutate(
        'df' = .data$nu - P + 1,
        'Sigma' = .data$Sigma / (.data$lambda * .data$df) ) %>%
      dplyr::group_by(.data$Draw, .data$Group) %>%
      dplyr::reframe(
        mvtnorm::rmvt(
          n = nb_sample,
          sigma = matrix(.data$Sigma, ncol = P, nrow = P),
          df = unique(.data$df),
          delta = unique(.data$mu)
        ) %>%
          `colnames<-`(unique(.data$Peptide)) %>%
          tibble::as_tibble(),
        'ID_sample' = 1:nb_sample
      ) %>%
      tidyr::pivot_longer(!c('Draw', 'Group', 'ID_sample'),
                          names_to = 'Peptide',
                          values_to = 'Sample') %>%
      dplyr::group_by(.data$Peptide, .data$Group, .data$ID_sample) %>%
      dplyr::reframe('Sample' = mean(.data$Sample)) %>%
      dplyr::select(- .data$ID_sample)
  }

  if('beta' %in% names(posterior)){
    dist = posterior %>%
      dplyr::mutate('sigma' = sqrt(.data$beta/(.data$lambda * .data$alpha))) %>%
      tidyr::uncount(nb_sample) %>%
      dplyr::reframe(
        .data$Peptide,
        .data$Group,
        'Sample' = extraDistr::rlst(
          n = dplyr::n(),
          df = 2 * .data$alpha,
          mu = .data$mu,
          sigma = .data$sigma)
      )
  }

  return(dist)
}
