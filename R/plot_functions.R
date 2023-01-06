#' Plot the difference of means' posterior distribution
#'
#' Display the posterior distribution of the difference of means between two
#' groups for a specific peptide. If only one group is provide, the function
#' display the posterior distribution of the mean for this specific group
#' instead. The function provides additional tools to represent information to
#' help inference regarding the difference between groups (reference at 0 on the
#' x-axis, probability of \code{group1} > \code{group2} and conversely).
#'
#' @param post_distrib A data frame, typically coming from a post_mean_*()
#'    function, containing the following columns: \code{Peptide}, \code{Group}
#'    and \code{Mean}. This argument should contain the empirical posterior
#'    distributions to be displayed.
#' @param group1 A character string, corresponding to the name of the group
#'    for which we plot the posterior distribution of the mean. If \code{group2}
#'    is provided, the posterior difference of the groups is displayed instead.
#' @param group2 A character string, corresponding to the name of the group
#'    we want to compare to \code{group1}. If NULL (default), only the posterior
#'    distribution of the mean for \code{group1} is displayed.
#' @param peptide A character string, corresponding to the name of the peptide
#'    for which we plot the posterior distribution of the mean. If NULL
#'    (default), all peptides contained in \code{post_distrib} are aggregated.
#' @param prob_CI A number, between 0 and 1, corresponding the level of the
#'    Credible Interval (CI), represented as side regions (in red) of the
#'    posterior distribution. The default value (0.95) display the 95% CI,
#'    meaning that the central region (in blue) contains 95% of the probability
#'    distribution of the mean.
#' @param show_prob A boolean, indicating whether we display the label of
#'    probability comparisons between the two groups.
#' @param mean_bar A boolean, indicating whether we display the vertical bar
#'    corresponding to 0 on the x-axis (when comparing two groups), of the mean
#'    value of the distribution (when displaying a unique group).
#' @param index_group1 A character string, used as the index of \code{group1} in
#'    the legends. If NULL (default), \code{group1} is used.
#' @param index_group2 A character string, used as the index of \code{group2} in
#'    the legends. If NULL (default), \code{group2} is used.
#'
#' @return Plot of the required posterior distribution.
#' @export
#'
#' @examples
#' TRUE
plot_posterior = function(
    post_distrib,
    group1,
    group2 = NULL,
    peptide = NULL,
    prob_CI = 0.95,
    show_prob = TRUE,
    mean_bar = TRUE,
    index_group1 = NULL,
    index_group2 = NULL
    ){
    ## Retrieve the name of peptides in the 'post_distrib' argument if needed
    if(peptide %>% is.null()){
      peptide = post_distrib$Peptide %>% unique()
    }

    ## If we have multiple values in 'peptide', warn the user
    if(length(peptide) > 1){
      cat("The 'peptide' argument contains multiple values. This might cause",
          "wrong results or unexpected graphs.")
    }

    ## Extract the distribution of group1
    db = post_distrib %>%
      dplyr::filter(.data$Peptide %in% peptide) %>%
      dplyr::filter(.data$Group == group1) %>%
      dplyr::pull(Mean)

    ## Compte the mean of the distribution to display as a vertical bar
    bar = mean(db)

    ## If group2 is provided, display the difference of posterior distributions
    if(!is.null(group2)){
      ## Extract the distribution of group1
      db2 = post_distrib %>%
        dplyr::filter(.data$Peptide %in% peptide) %>%
        dplyr::filter(.data$Group == group2) %>%
        dplyr::pull(.data$Mean)

      ## Redefine reference distribution as the difference of group1 and group2
      db = db - db2
      bar = 0
      }

  ## Create a density from the dataset
  dens <- density(db, n = 5000)

  ## Compute the Credible Interval with the appropriate probability level
  CI = quantile(db, prob= c( (1 - prob_CI)/2, (1 + prob_CI)/2 ) )

  ## Format the dataset for the subsequent plot
  db_plot = tibble::tibble(x = dens$x, y = dens$y) %>%
    dplyr::mutate(quant = factor(findInterval(x, CI)))

  ## Define the name of index for the label of group1
  if(index_group1 %>% is.null()){
    index_group1 = group1
  }
  ## Define the name of index for the label of group2
  if(index_group2 %>% is.null()){
    index_group2 = group2
  }

  gg = ggplot2::ggplot(db_plot) +
    ggplot2::geom_ribbon(
      ggplot2::aes(x = x, y = y, ymin=0, ymax=y, fill = quant)
      ) +
    ggplot2::ylab('Density') +
    ggplot2::scale_fill_manual(values=c("#F8B9C5", "#AFC0E3", "#F8B9C5")) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position="none")

  if(mean_bar){
    gg = gg + ggplot2::geom_vline(xintercept = bar, color = 'red')
  }

  if(group2 %>% is.null()){
    ## Add the adequate label for the x-axis
    gg = gg + ggplot2::xlab( bquote(mu[.(index_group1)]) )
  } else {
      ## Add the adequate label for the x-axis
      gg = gg + ggplot2::xlab(bquote(mu[.(index_group1)] - mu[.(index_group2)]))
      ## Add probabilities of the group comparison if required
     if( show_prob == TRUE ){
      p_inf = (sum(db<0)/length(db)) %>% round(2) %>% as.character()
      p_sup = (sum(db>0)/length(db)) %>% round(2) %>% as.character()
      exp_l = bquote(P(mu[.(index_group1)] <= mu[.(index_group2)]) == .(p_inf))
      exp_r = bquote(P(mu[.(index_group1)] >= mu[.(index_group2)]) == .(p_sup))

      gg = gg +
        ggplot2::geom_label(
          data = tibble::tibble(bar = bar),
          ggplot2::aes(
           x = bar,
            y = Inf,
            label = deparse(exp_l)
            ),
          parse = T,
          size = 4,
          hjust=1,
          vjust=1) +
          ggplot2::geom_label(
            data = tibble::tibble(bar = bar),
            ggplot2::aes(
              x = bar,
              y = Inf,
              label = deparse(exp_r)
            ),
            parse = T,
            size = 4,
            hjust=0,
            vjust=1)
      }
    }
  return(gg)
}
