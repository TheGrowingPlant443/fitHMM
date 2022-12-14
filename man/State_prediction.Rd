% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/State_prediction.R
\name{State_prediction}
\alias{State_prediction}
\title{State Prediction Function}
\usage{
State_prediction(HMM, observations, steps = 1)
}
\arguments{
\item{HMM}{An object of class \code{HMM} or \code{xHMM}.}

\item{observations}{A numeric vector containing the observed sequence.}

\item{steps}{Numeric of length 1, the number of steps.}
}
\value{
Vector of the probabilities of being in each hidden state \code{steps}-time
units in the future.
}
\description{
Calculates the probability of being in each hidden state at time 
T + steps, where T is the length of the observed sequence, i.e. 
the future hidden states probability.
}
