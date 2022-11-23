#' @import expm
#' @title State Prediction Function
#' @description Calculates the probability of being in each hidden state at time 
#' T + steps, where T is the length of the observed sequence, i.e. 
#' the future hidden states probability.
#' 
#' 
#' @usage State_prediction(HMM, observations, steps = 1)
#' 
#' @param HMM An object of class \code{HMM} or \code{xHMM}.
#' @param observations A numeric vector containing the observed sequence.
#' @param steps Numeric of length 1, the number of steps.
#'
#' @return Vector of the probabilities of being in each hidden state \code{steps}-time
#' units in the future.
#' 
#' @export
State_prediction <- function(HMM, observations,steps=1){
  AlphaBetaLikelihood <- ForwardBackwardAlgorithm(HMM, observations) 
  Alpha <- AlphaBetaLikelihood$Alpha
  Beta <- AlphaBetaLikelihood$Beta
  Likelihood <- AlphaBetaLikelihood$Likelihood_forw
  State_predict_after_T = (Alpha[,dim(Alpha)[2]] %*% (expm::`%^%`(HMM$TPM,steps)))/Likelihood
  return(State_predict_after_T)
}
