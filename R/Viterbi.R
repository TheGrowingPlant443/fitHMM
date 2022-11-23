### Generic Function for Viterbi Decoding --------------------------------------

#' @title Viterbi Decoding for HMM Objects
#'
#' @description Generic function for performing Viterbi, i.e. global, decoding of
#' an HMM and an observed sequence from said HMM.
#' 
#' @usage Viberbi(HMM_obj, ...)
#' 
#' ## S3 methods for PoissonHMM, BernoulliHMM, NormalHMM, xHMM:
#' Viterbi(HMM_obj, obs, ...)
#' 
#'
#' @param HMM_obj an object of class \code{HMM} or \code{xHMM}.
#' @param obs a numeric vector containing the observed sequence.
#'
#' @details Viterbi decoding finds the most probable sequence of hidden states
#' given the observed sequence \code{obs}.
#' The function uses scaling to avoid underflow.
#' 
#' @return Returns a numeric vector containing the predicted hidden states that
#' emitted the observed sequence \code{obs}. The hidden states are numbered
#' according the order of the emission parameters \code{HMM_obj$emis_param}.
#' 
#' @export
Viterbi <- function(HMM_obj, ...){
  UseMethod("Viterbi")
}

### Viterbi Decoding for Poisson -----------------------------------------------
#' @export
Viterbi.PoissonHMM <- function(HMM_obj, obs, ...){
  is.numeric(obs) || stop("obs must be a numeric vector")
  TT <- length(obs) # Length of observed sequence
  m <- HMM_obj$nstates # nr. of hidden states
  lambda <- HMM_obj$emis_param$lambda # Rates of hidden states
  Gam <- HMM_obj$TPM # TPM of hidden Markov chain
  ## Create matrices:
  xi <- matrix(0, nrow = m, ncol = TT)
  A <- matrix(0, nrow = m, ncol = TT-1) # Arrow matrix
  ## Initialize:
  element <- HMM_obj$initial_probs*dpois(obs[1], lambda)
  xi[,1] <- element/sum(element)
  ## Recursion:
  for (t in 2:TT){
    element <- dpois(obs[t], lambda)*apply(xi[,t-1]*Gam, 2, max)
    xi[,t] <- element/sum(element)
    A[,t-1] <- apply(xi[,t-1]*Gam, 2, which.max)
  }
  ## Backtrack:
  decode <- numeric(TT)
  decode[TT] <- which.max(xi[,TT]) # Final element of most probable sequence
  for (t in (TT-1):1){
    decode[t] <- A[decode[t+1], t] # Element that lead to next
  }
  return(decode)
}

### Viterbi Decoding for Bernoulli ---------------------------------------------
#' @export
Viterbi.BernoulliHHM <- function(HMM_obj, obs, ...){
  is.numeric(obs) || stop("obs must be a numeric vector")
  TT <- length(obs) # Length of observed sequence
  m <- HMM_obj$nstates # nr. of hidden states
  prob <- HMM_obj$emis_param$prob # pbbs of hidden states
  Gam <- HMM_obj$TPM # TPM of hidden Markov chain
  ## Create matrices:
  xi <- matrix(0, nrow = m, ncol = TT)
  A <- matrix(0, nrow = m, ncol = TT-1) # Arrow matrix
  ## Initialize:
  element <- HMM_obj$initial_probs*dbinom(obs[1], size = 1, prob = prob)
  xi[,1] <- element/sum(element)
  ## Recursion:
  for (t in 2:TT){
    element <- dbinom(obs[t], size = 1, prob = prob)*apply(xi[,t-1]*Gam, 2, max)
    xi[,t] <- element/sum(element)
    A[,t-1] <- apply(xi[,t-1]*Gam, 2, which.max)
  }
  ## Backtrack:
  decode <- numeric(TT)
  decode[TT] <- which.max(xi[,TT]) # Final element of most probable sequence
  for (t in (TT-1):1){
    decode[t] <- A[decode[t+1], t] # Element that lead to next
  }
  return(decode)
}


### Viterbi Decoding for Normal ------------------------------------------------
#' @export
Viterbi.NormalHMM <- function(HMM_obj, obs, ...){
  is.numeric(obs) || stop("obs must be a numeric vector")
  TT <- length(obs) # Length of observed sequence
  m <- HMM_obj$nstates # nr. of hidden states
  mu <- HMM_obj$emis_param$mean
  sig <- HMM_obj$emis_param$sd
  Gam <- HMM_obj$TPM # TPM of hidden Markov chain
  ## Create matrices:
  xi <- matrix(0, nrow = m, ncol = TT)
  A <- matrix(0, nrow = m, ncol = TT-1) # Arrow matrix
  ## Initialize:
  element <- HMM_obj$initial_probs*dnorm(obs[1], mean = mu, sd = sig)
  xi[,1] <- element/sum(element)
  ## Recursion:
  for (t in 2:TT){
    element <- dnorm(obs[t], mean = mu, sd = sig)*apply(xi[,t-1]*Gam, 2, max)
    xi[,t] <- element/sum(element)
    A[,t-1] <- apply(xi[,t-1]*Gam, 2, which.max)
  }
  ## Backtrack:
  decode <- numeric(TT)
  decode[TT] <- which.max(xi[,TT]) # Final element of most probable sequence
  for (t in (TT-1):1){
    decode[t] <- A[decode[t+1], t] # Element that lead to next
  }
  return(decode)
}


### Viterbi decoding for xHMM --------------------------------------------------
#' @export
Viterbi.xHMM <- function(HMM_obj, obs, ...){
  is.numeric(obs) || stop("obs must be a numeric vector")
  TT <- length(obs) # Length of observed sequence
  m <- HMM_obj$nstates # nr. of hidden states
  param <- HMM_obj$emis_param # Emission parameters
  Gam <- HMM_obj$TPM # TPM of hidden Markov chain
  dist <- HMM_obj$distribution
  ## Create matrices:
  xi <- matrix(0, nrow = m, ncol = TT)
  A <- matrix(0, nrow = m, ncol = TT-1) # Arrow matrix
  ## Initialize:
  element <- HMM_obj$initial_probs*do.call(dist, append(list(obs[1]), param))
  xi[,1] <- element/sum(element)
  ## Recursion:
  for (t in 2:TT){
    element <- do.call(dist, append(list(obs[t]), param))*apply(xi[,t-1]*Gam, 2, max)
    xi[,t] <- element/sum(element)
    A[,t-1] <- apply(xi[,t-1]*Gam, 2, which.max)
  }
  ## Backtrack:
  decode <- numeric(TT)
  decode[TT] <- which.max(xi[,TT]) # Final element of most probable sequence
  for (t in (TT-1):1){
    decode[t] <- A[decode[t+1], t] # Element that lead to next
  }
  return(decode)
}

