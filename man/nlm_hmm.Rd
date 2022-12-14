% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/general_nlm.R
\name{nlm_hmm}
\alias{nlm_hmm}
\title{Numeric Parameter Estimation of HMM}
\usage{
nlm_hmm(HMM, observations)
}
\arguments{
\item{HMM}{object of class \code{xHMM} containing the initial parameters.}

\item{observations}{numeric vector containing the observed sequence.}
}
\value{
Returns a new object of class \code{xHMM} that is fitted to the
observations according to non-linear minimization approach.
}
\description{
Uses \code{nlm} to minimize the negative log likelihood of an HMM,
which corresponds to find the parameters maximizing the likelihood.
}
