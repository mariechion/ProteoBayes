% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_functions.R
\name{plot_multi_diff}
\alias{plot_multi_diff}
\title{Plot multivariate comparison between groups}
\usage{
plot_multi_diff(proba_diff, cumulative = FALSE)
}
\arguments{
\item{proba_diff}{A tibble, typically coming from the
\code{multi_identify_diff} function, containing probability distribution
of the number of larger peptides in each one-to-one group comparisons.}
}
\value{
A graph (or a matrix of graphs) displaying the multivariate
differential inference summary between groups
}
\description{
Graphical representation of inference in a multivariate difference
analysis context. The plotted empirical distribution represents, for each
one-to-one group comparison, the probability of the number of elements for
which the mean of a peptide is larger in a given group. This provides visual
evidence on whether two groups are differential or not, with adequate
uncertainty quantification.
}
\examples{
TRUE
}
