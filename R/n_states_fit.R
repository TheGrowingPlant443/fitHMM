#' @title Fit n-State HMM
#' 
#' @description Fit an HMM with a standard emission distribution and a specified
#' number of hidden states.
#' 
#' @usage Nstate_HMM(distribution = "Normal", data, nstates = 2)
#' 
#' @param distribution a character of length one specifying the emission
#' distribution, must be one of the following: "Poisson", "Bernoulli", "Normal".
#' @param data a numeric vector containing the observed data.
#' @param nstates a numeric of length one, the size of the hidden state space.
#' 
#' @details \code{Nstate_HMM} calls \code{HMM_EM} with initial parameters chosen
#' using summary statistics of \code{data}.
#'
#' @return Returns an object of class \code{HMM} containing parameter estimates
#' for the fitted model.
#' 
#' @seealso \code{\link{HMM_EM}} for the EM-algorithm.
#' 
#' @example
#' ## Fit 3-state Poisson HMM to Poisson2State:
#' nstate_HMM(distribution = "Poisson", fitHMM::Poisson2State, nstates = 3)
#' 
#' 
#' @export
Nstate_HMM <- function(distribution="Normal", data, nstates=2){
  Time = length(data)
  nstates <= Time || stop("number of states must be less than number of datapoints")
  means = quantile(data, (1:nstates-1/2)*1/nstates)
  if (distribution=="Normal"){
    parameters = list(mean=means, sd=rep(1, nstates))
  }
  else if (distribution=="Poisson"){
    parameters = list(lambda=means)
  }
  else if (distribution=="Bernoulli"){
    parameters = list(prob=(1:nstates-1/2)*1/nstates)
  }
  else{
    stop("distribution needs to be one of: Normal, Poisson, Bernoulli")
  }
  model = HMM(distribution = distribution,
              emis_param = parameters,
              initial_probs = rep(1/nstates, nstates),
              TPM = matrix(1/nstates, nstates, nstates))
  model = HMM_EM(model, data)
  return(model)
}