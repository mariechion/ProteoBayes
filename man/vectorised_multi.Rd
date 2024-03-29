% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_posteriors.R
\name{vectorised_multi}
\alias{vectorised_multi}
\title{Vectorised version of multi_posterior_mean()}
\usage{
vectorised_multi(data, mu_0 = NULL, lambda_0 = 1, Sigma_0 = NULL, nu_0 = 10)
}
\arguments{
\item{data}{A tibble or data frame containing imputed data sets for all
groups. Required columns: \code{Peptide}, \code{Group}, \code{Sample},
\code{Output}. If missing data have been estimated from multiple
imputations, each imputation should be identified in an optional
\code{Draw} column.}

\item{mu_0}{A vector, corresponding to the prior mean. If NULL, all groups
are initialised with the same empirical mean for each peptide.}

\item{lambda_0}{A number, corresponding to the prior covariance scaling
parameter.}

\item{Sigma_0}{A matrix, corresponding to the prior covariance parameter.
If NULL, the identity matrix will be used by default.}

\item{nu_0}{A number, corresponding to the prior degrees of freedom.}
}
\value{
A tibble providing the parameters of the posterior t-distribution for
the mean of the considered groups for each peptide.
}
\description{
Alternative vectorised version, highly efficient when nb_peptide < 30.
}
\examples{
TRUE
}
