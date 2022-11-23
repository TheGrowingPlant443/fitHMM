### Documentation of Nbinom2State data set

#' @title Nbinom2State Dataset
#' 
#' @description An observation of length 100 simulated from a two-state negative binomial HMM.
#' 
#' @details The data is simulated with the parameters:
#' \describe{
#' \item{sizes:}{10, 100}
#' \item{probs:}{0.25, 0.30}
#' \item{initial probability for state one:}{0.5}
#' \item{transition probability from one to one:}{0.6}
#' \item{transition probability from two to two:}{0.8}
#' }
#' 
#' @docType data
#' 
#' @format a numeric vector of length 100.
"Nbinom2State"