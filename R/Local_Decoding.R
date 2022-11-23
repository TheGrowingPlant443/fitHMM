#' @import expm

#' @title Local Decoding
#' 
#' @description Calculate the probability of being in each state to time
#' time T.
#' 
#' @usage Local_decoding(HMM, observations)
#'
#' @param HMM Your specified hidden Markov model by using the function \code{HMM}
#' or \code{xHMM}
#' @param observations Your Data
#' 
#'
#' @return Vector of probability of being in each state
#' @export Local_decoding
Local_decoding <- function(HMM, observations){
  AlphaBetaLikelihood <- ForwardBackwardAlgorithm(HMM, observations) 
  Alpha <- AlphaBetaLikelihood$Alpha
  Beta <- AlphaBetaLikelihood$Beta
  Likelihood <- AlphaBetaLikelihood$Likelihood_forw
  Local_decoding_matrix <- matrix(0,nrow = dim(Alpha)[1], ncol = dim(Alpha)[2])
  for (i in 1:dim(Alpha)[1]){
    for (j in 1:dim(Alpha)[2]) {
      Local_decoding_matrix[i,j] <- (Alpha[i,j]*Beta[i,j])/Likelihood
    }
  }
  Local_decoding_index <- c()
  Local_decoding_probability <- c()
  for (j in 1:dim(Alpha)[2]) {
    Local_decoding_index <- c(Local_decoding_index,
                              which.max(Local_decoding_matrix[,j]))
    Local_decoding_probability <- c(Local_decoding_probability,
                                    max(Local_decoding_matrix[,j]))
  }
  return(list(Local_decoding_matrix=Local_decoding_matrix,
              Local_decoding_index = Local_decoding_index, 
              Local_decoding_probability = Local_decoding_probability
  ))
}

#' @title Local_decoding_with_scaling
#'
#' @description Calculate the probability of being in each state to time
#' time T with scaling.
#' 
#' @usage Local_decoding_with_scaling(HMM, observations)
#' 
#' @param HMM Your specified hidden Markov model by using the function \code{HMM}
#' or \code{xHMM}
#' @param observations Your Data
#' 
#'
#' @return Vector of probability of being in each state with scaling
#' 
#' @export Local_decoding_with_scaling
Local_decoding_with_scaling <- function(HMM, observations){
  AlphaBetaLikelihood <- ForwardBackwardAlgorithm(HMM, observations) 
  Alpha <- AlphaBetaLikelihood$Alpha_scaled
  Beta <- AlphaBetaLikelihood$Beta_scaled
  Likelihood <- AlphaBetaLikelihood$LogLikelihood_forw
  Local_decoding_matrix <- matrix(0,nrow = dim(Alpha)[1], ncol = dim(Alpha)[2])
  for (i in 1:dim(Alpha)[1]){
    for (j in 1:dim(Alpha)[2]) {
      Local_decoding_matrix[i,j] <- exp((Alpha[i,j]+Beta[i,j]-Likelihood)/1000)
    }
  }
  Local_decoding_index <- c()
  Local_decoding_probability <- c()
  for (j in 1:dim(Alpha)[2]) {
    Local_decoding_index <- c(Local_decoding_index,
                              which.max(Local_decoding_matrix[,j]))
    Local_decoding_probability <- c(Local_decoding_probability,
                                    max(Local_decoding_matrix[,j]))
  }
  return(list(Local_decoding_matrix=Local_decoding_matrix,
              Local_decoding_index = Local_decoding_index, 
              Local_decoding_probability = Local_decoding_probability
              ))
}

