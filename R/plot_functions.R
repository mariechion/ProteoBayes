#' Plot the posterior distribution of the difference of means
#'
#' Display the posterior distribution of the difference of means between two
#' groups for a specific peptide. If only one group is provide, the function
#' display the posterior distribution of the mean for this specific group
#' instead. The function provides additional tools to represent information to
#' help inference regarding the difference between groups (reference at 0 on the
#' x-axis, probability of \code{group1} > \code{group2} and conversely).
#'
#' @param sample_distrib A data frame, typically coming from the
#'    \code{sample_distrib()} function, containing the following columns:
#'    \code{Peptide}, \code{Group} and \code{Sample}. This argument should
#'    contain the empirical posterior distributions to be displayed.
#' @param group1 A character string, corresponding to the name of the group
#'    for which we plot the posterior distribution of the mean. If NULL
#'    (default), the first group appearing in \code{sample_distrib} is
#'    displayed. If \code{group2} is provided, the posterior difference of the
#'    groups is displayed instead.
#' @param group2 A character string, corresponding to the name of the group
#'    we want to compare to \code{group1}. If NULL (default), only the posterior
#'    distribution of the mean for \code{group1} is displayed.
#' @param peptide A character string, corresponding to the name of the peptide
#'    for which we plot the posterior distribution of the mean. If NULL
#'    (default), only the first appearing in \code{sample_distrib} is displayed.
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
plot_distrib = function(
    sample_distrib,
    group1 = NULL,
    group2 = NULL,
    peptide = NULL,
    prob_CI = 0.95,
    show_prob = TRUE,
    mean_bar = TRUE,
    index_group1 = NULL,
    index_group2 = NULL
    ){

  ## Retrieve the name of the first group in 'sample_distrib' if needed
  if(group1 %>% is.null()){
    group1 = sample_distrib$Group[1]
  }

  ## Retrieve the name of peptides in the 'sample_distrib' argument if needed
  if(peptide %>% is.null()){
    peptide = sample_distrib$Peptide %>% unique()
  }

  ## If we have multiple values in 'peptide', warn the user
  if(length(peptide) > 1){
    peptide = peptide[1]
  }

  ## Extract the distribution of group1
  db = sample_distrib %>%
    dplyr::filter(.data$Peptide %in% peptide) %>%
    dplyr::filter(.data$Group == group1) %>%
    dplyr::pull(.data$Sample)

  ## Compute the mean of the distribution to display as a vertical bar
  bar = mean(db)

  ## If group2 is provided, display the difference of posterior distributions
  if(!is.null(group2)){
    ## Extract the distribution of group1
    db2 = sample_distrib %>%
      dplyr::filter(.data$Peptide %in% peptide) %>%
      dplyr::filter(.data$Group == group2) %>%
      dplyr::pull(.data$Sample)

    ## Redefine reference distribution as the difference of group1 and group2
    db = db - db2
    bar = 0
  }

  ## Create a density from the dataset
  dens <- stats::density(db, n = 5000)

  ## Compute the Credible Interval with the appropriate probability level
  CI = stats::quantile(db, prob= c( (1 - prob_CI)/2, (1 + prob_CI)/2 ) )

  ## Format the dataset for the subsequent plot
  db_plot = tibble::tibble('x' = dens$x, 'y' = dens$y) %>%
    dplyr::mutate('quant' = factor(findInterval(.data$x, CI)))

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
      ggplot2::aes(x = .data$x,
                   y = .data$y,
                   ymin=0,
                   ymax=.data$y,
                   fill = .data$quant)
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

#' Plot multivariate comparison between groups
#'
#' Graphical representation of inference in a multivariate difference
#' analysis context. The plotted empirical distribution represents, for each
#' one-to-one group comparison, the probability of the number of elements for
#' which the mean of a peptide is larger in a given group. This provides visual
#' evidence on whether two groups are differential or not, with adequate
#' uncertainty quantification.
#'
#' @param multi_diff A tibble, typically coming from the
#'    \code{multi_identify_diff} function, containing probability distribution
#'    (or the cumulative distribution) of the number of larger Peptides in each
#'     one-to-one group comparisons.
#'
#' @param cumulative A boolean, indicating whether the cumulative distribution
#'    or the original probability distribution should be displayed.
#' @param plot_mean A boolean, indicating whether the graph for difference of
#'    means between all groups for each Peptide should be displayed.
#'
#' @returns A graph (or a matrix of graphs) displaying the multivariate
#'      differential inference summary between groups
#' @export
#'
#' @examples
#' TRUE
plot_multi_diff = function(
  multi_diff,
  plot_mean = TRUE,
  cumulative = FALSE
  ){

  proba_diff = multi_diff$Diff_proba
  mean_diff = multi_diff$Diff_mean

  ## Initialise the list of graphs to be plotted
  gg = list()

  ## Initialise a counter for the number of graphs
  counter = 0

  ## Get the list of Groups
  list_groups = proba_diff$Group1 %>%
    union(proba_diff$Group2) %>%
    unique()

  ## Initialise the layout matrix for displaying multiple graphs
  layout_matrix = matrix(NA,
                         nrow=length(list_groups)-1,
                         ncol=length(list_groups)-1
                         )

  ## Initialise the list of remaining groups to avoid duplicated comparisons
  list_remaining_groups = list_groups

  ## Loop over all one-by-one group comparisons
  for(i in list_groups){

    ## Remove current group1 from the list of group2
    list_remaining_groups = list_remaining_groups[-1]

    for(j in list_remaining_groups){
      db_plot = proba_diff %>%
        dplyr::filter(.data$Group1 == i, .data$Group2 == j)

      ## Increment the counter and position it in the layout matrix
      counter = counter + 1
      layout_matrix[i,j-1] = counter

      ## Get the index of Group1 for the legend
      index_group1 = db_plot %>%
        dplyr::pull(.data$Group1) %>%
        unique() %>%
        as.character()

      ## Get the index of Group2 for the legend
      index_group2 = db_plot %>%
        dplyr::pull(.data$Group2) %>%
        unique() %>%
        as.character()

      ## Get the list of number of peptides for the x-axis
      nb_peptides = db_plot$Nb_peptides

      ## Display whether the probability distribution or its cumulative
      if('Proba' %in% names(db_plot)){
        gg[[counter]] = ggplot2::ggplot(db_plot) +
          ggplot2::geom_bar(
            ggplot2::aes(x = .data$Nb_peptides, y = .data$Proba),
            stat = 'identity',
            fill = '#AFC0E3'
          ) +
          ggplot2::geom_vline(ggplot2::aes(xintercept = max(nb_peptides)/2),
                              color = 'red',
                              linetype = 'dashed') +
          ggplot2::theme_classic() +
          ggplot2::xlab(
            bquote( paste('Number of peptides i where ',
                          mu[.(index_group1)]^i > mu[.(index_group2)]^i))) +
          ggplot2::ylab(bquote(Probability)) +
          ggplot2::xlim(c(min(nb_peptides)-0.5, max(nb_peptides)+0.5))

      } else {
        gg[[counter]] = ggplot2::ggplot(db_plot) +
          ggplot2::geom_line(
            ggplot2::aes(x = .data$Nb_peptides, y = .data$Cumul_proba),
            col = '#AFC0E3') +
          ggplot2::theme_classic() +
          ggplot2::xlab(
            bquote(paste('Number of peptides i where ',
                         mu[.(index_group1)]^i > mu[.(index_group2)]^i))) +
          ggplot2::ylab(bquote('Cumulative probability')) +
          ggplot2::xlim(c(min(nb_peptides), max(nb_peptides)))
      }
      ## Add overlapping coefficient in title if provided in 'multi_diff'
      if('Overlap_coef' %in% names(multi_diff)){

        coef = multi_diff$Overlap_coef %>%
          dplyr::filter(.data$Group1 == i, .data$Group2 == j) %>%
          dplyr::pull(.data$Overlap_coef)

        gg[[counter]] = gg[[counter]] +
          ggtitle(bquote(paste('Overlapping coefficient:', .(coef))))
      }
    }
  }

  if(plot_mean){

    if(length(list_groups) == 2){
     layout_matrix = as.matrix(c(1,2))
    } else{
     layout_matrix[length(list_groups)-1, 1] = counter+1
    }

    gg[[counter+1]] = ggplot2::ggplot(mean_diff) +
      ggplot2::geom_point(ggplot2::aes(x = .data$Mean,
                                       y = .data$Peptide,
                                       col = factor(.data$Group))) +
      ggplot2::xlab(bquote('Posterior Mean')) +
      ggplot2::theme_classic() +
      ggplot2::scale_colour_discrete(name = "Group")
  }

  gridExtra::grid.arrange(grobs=gg, layout_matrix = layout_matrix) %>%
    return()
}
