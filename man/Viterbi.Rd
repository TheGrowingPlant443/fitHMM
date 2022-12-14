% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Viterbi.R
\name{Viterbi}
\alias{Viterbi}
\title{Viterbi Decoding for HMM Objects}
\usage{
Viberbi(HMM_obj, ...)

## S3 methods for PoissonHMM, BernoulliHMM, NormalHMM, xHMM:
Viterbi(HMM_obj, obs, ...)
}
\arguments{
\item{HMM_obj}{an object of class \code{HMM} or \code{xHMM}.}

\item{obs}{a numeric vector containing the observed sequence.}
}
\value{
Returns a numeric vector containing the predicted hidden states that
emitted the observed sequence \code{obs}. The hidden states are numbered
according the order of the emission parameters \code{HMM_obj$emis_param}.
}
\description{
Generic function for performing Viterbi, i.e. global, decoding of
an HMM and an observed sequence from said HMM.
}
\details{
Viterbi decoding finds the most probable sequence of hidden states
given the observed sequence \code{obs}.
The function uses scaling to avoid underflow.
}
