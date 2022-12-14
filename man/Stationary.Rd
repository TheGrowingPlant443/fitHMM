% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Model Checking.R
\name{Stationary}
\alias{Stationary}
\title{Find Stationary Distribution of HMM}
\usage{
Stationary(HMM)
}
\arguments{
\item{HMM}{the \code{HMM} or \code{xHMM} object to determine the stationary
distribution of.}
}
\value{
a \code{HMM} or \code{xHMM} object identical to the input except the
\code{initial_probs} parameter has been changed to the stationary distribution
of the (hidden) Markov chain.
}
\description{
Find Stationary Distribution of HMM
}
