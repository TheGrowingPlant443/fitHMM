# Forward and Backward Algorithm ----
#' Forward and Backward Algorithm for HMMs
#' 
#' @description Runs the forward and backward algorithm for an HMM and observation
#' to calculate the forward and backward probabilities.
#' 
#' @usage ForwardBackwardAlgorithm(HMM, observations)
#'
#' @param HMM object of class \code{HMM} or \code{xHMM} containing the parameters
#' of the HMM.
#' @param observations a numeric vector of observations.
#' 
#' @details The function calculates the forward probabilities,
#' \deqn{\alpha_t(i) = P(c_t = i, x_1,\dots,x_t),}
#' and the backward probabilities,
#' \deqn{\beta_t(i) = P(x_{t+1},\dots,x_T\vert c_t=i) \quad t=1,\dots,T-1,}
#' where \eqn{x_1,\dots,x_T} denotes the observed sequence and \eqn{c_1,\dots,c_T}
#' denotes the hidden Markov chain.
#' 
#' The forward or the backward probabilities can be used to calculate the likelihood
#' function and the state probabilities, i.e. \eqn{P(c_t=i\vert x)}.
#' 
#' @return Returns a list with elements:
#' \describe{
#' \item{Alpha}{a matrix containing the forward probabilities with time along the
#' columns and hidden states along the rows.}
#' \item{Alpha_scaled}{contains a scaled version of the forward probabilities.}
#' \item{Beta}{a matrix containing the backward probabilities with time along the
#' columns and hidden states along the rows.}
#' \item{Beta_scaled}{contains a scaled version of the backward probabilities.}
#' \item{Likelihood_forw}{a numeric containing the likelihood function evaluated
#' in \code{observation} and parameters from \code{HMM}. Computed using the forward
#' probabilities.}
#' \item{Likelihood_back}{a numeric containing the likelihood function evaluated
#' in \code{observation} and parameters from \code{HMM}. Computed using the backward
#' probabilities. }
#' \item{P_mat}{a matrix containing the state probabilities with time along the columns
#' and hidden states along the rows.}
#' }
#' 
#' 
#' @export
ForwardBackwardAlgorithm <- function(HMM, observations){
  Time = length(observations)
  if (class(HMM)[1]=="xHMM"){
    ddist = HMM$distribution
  }
  else{
    ddist = list(Poisson=dpois, Normal=dnorm, Bernoulli=function(x, prob){dbinom(x, 1, prob)})[HMM$distribution][[1]]
  }
  P_mat = matrix(nrow=HMM$nstates, ncol=Time)
  for (i in 1:Time){
    P_mat[,i] = do.call(ddist, append(list(observations[i]), HMM$emis_param))
  }
  Alpha = matrix(nrow=HMM$nstates, ncol=Time)
  Beta = matrix(nrow=HMM$nstates, ncol=Time)
  Alpha[,1] = matrix(HMM$initial_probs * P_mat[,1])
  Beta[,Time] = rep(1, HMM$nstates)
  Alpha_scaled = Alpha
  Beta_scaled = Beta
  Beta_temp = Beta
  AS = c()
  AS[1] = sum(Alpha[,1])
  BS = log(length(HMM$nstates))
  for (t in 2:Time){
    t_1 = Time-t+1
    Alpha[,t] = Alpha[,t-1] %*% HMM$TPM * P_mat[,t]
    Alpha_scaled[,t] = (Alpha_scaled[,t-1] %*% HMM$TPM * P_mat[,t]) / AS[t-1]
    Beta[,t_1] = HMM$TPM %*% (Beta[,t_1+1] * P_mat[,t_1+1])
    Beta_temp[,t_1] = HMM$TPM %*% (Beta_temp[,t_1+1] * P_mat[,t_1+1])
    Beta_scaled[,t_1] = log(Beta_temp[,t_1]) + BS
    Beta_temp[,t_1] = Beta_temp[,t_1]/sum(Beta_temp[,t_1])
    BS = BS + log(sum(Beta_temp[,t_1]))
    AS[t] = sum(Alpha_scaled[,t])
  }
  return(list(Alpha = Alpha, 
              Alpha_scaled = Alpha_scaled, 
              Beta = Beta,
              Beta_scaled = Beta_scaled,
              Likelihood_forw = exp(sum(log(AS))),
              LogLikelihood_forw = sum(log(AS)),
              Likelihood_back = sum(Beta[,1] * P_mat[,1] * HMM$initial_probs),
              LogLikelihood_back = sum(log(AS)),
              P_mat = P_mat))
}

# EM Algorithm - Poisson ----

#' @title Expectation Maximization Algorithm for HMM Objects
#' 
#' @description The EM-algorithm for objects of class \code{HMM}.
#' 
#' @usage HMM_EM(HMM, ...)
#' 
#' ## S3 methods for PoissonHMM, BernoulliHMM, NormalHMM:
#' HMM_EM(HMM, observations, ...)
#'
#' @param HMM an object of class \code{HMM} containing model and initial parameters.
#' @param observations a numeric vector containing the observed sequence.
#' 
#' @details The EM-algorithm, also called the Baum-Welch algorithm in the HMM context,
#' can be used to find approximate maximum likelihood estimates for the parameters.
#' The desired model, i.e. number of hidden states and emission distribution, is
#' specified in \code{HMM}, which contains this specification along with initial
#' parameters for the EM-algorithm.
#' 
#' @return Returns an object of class \code{HMM} containing the estimated parameters.
#' 
#' @seealso \code{\link{Nstate_HMM}} for easy fitting of an HMM with a standard emission
#' distribution.
#' 
#' @export
HMM_EM <- function(HMM, ...){
  UseMethod("HMM_EM")
}

#' @export
HMM_EM.PoissonHMM <- function(HMM, observations){
  Time = length(observations)
  while(T){
    FB_alg = ForwardBackwardAlgorithm(HMM, observations)
    N_ij = matrix(nrow=HMM$nstates, ncol=HMM$nstates)
    for (i in 1:HMM$nstates){
      for (j in 1:HMM$nstates){
        N_ij[i,j] = with(FB_alg, 
                         sum(Alpha[i,-Time] * Beta[j,-1] * P_mat[j,-1] * HMM$TPM[i,j] / Likelihood_forw))
      }
    }
    N_ia = with(FB_alg,
                rowSums(Alpha * Beta) / Likelihood_forw)
    N_aa = with(FB_alg,
                as.numeric((Alpha * Beta) %*% observations) / Likelihood_forw)
    new_initial_probs = with(FB_alg,
                             (Alpha[,1] * Beta[,1])/Likelihood_forw)
    new_TPM = N_ij / rowSums(N_ij)
    new_lambda = N_aa / N_ia
    if (sum(abs(new_TPM - HMM$TPM)) <= 0.0001 &&
        sum(abs(new_initial_probs - HMM$initial_probs)) <= 0.0001 &&
        sum(abs(new_lambda - HMM$emis_param$lambda)) <= 0.0001){
      return(HMM)
    }
    HMM$TPM = new_TPM
    HMM$initial_probs = new_initial_probs
    HMM$emis_param$lambda = new_lambda
  }
}

#' @export
HMM_EM.NormalHMM <- function(HMM, observations){
  Time = length(observations)
  while(T){
    FB_alg = ForwardBackwardAlgorithm(HMM, observations)
    N_ij = matrix(nrow=HMM$nstates, ncol=HMM$nstates)
    for (i in 1:HMM$nstates){
      for (j in 1:HMM$nstates){
        N_ij[i,j] = with(FB_alg, 
                         sum(Alpha[i,-Time] * Beta[j,-1] * P_mat[j,-1] * HMM$TPM[i,j] / Likelihood_forw))
      }
    }
    N_ia = with(FB_alg,
                rowSums(Alpha * Beta) / Likelihood_forw)
    N_aa = with(FB_alg,
                as.numeric((Alpha * Beta) %*% observations) / Likelihood_forw)
    new_initial_probs = with(FB_alg,
                             (Alpha[,1] * Beta[,1])/Likelihood_forw)
    new_TPM = N_ij / rowSums(N_ij)
    new_mean = N_aa / N_ia
    new_sd = with(FB_alg, as.numeric((Alpha * Beta) %*% (observations - new_mean)^2) / Likelihood_forw) / N_ia
    if (sum(abs(new_TPM - HMM$TPM)) <= 0.0001 &&
        sum(abs(new_initial_probs - HMM$initial_probs)) <= 0.0001 &&
        sum(abs(new_mean - HMM$emis_param$mean)) <= 0.0001 &&
        sum(abs(new_sd - HMM$emis_param$sd))){
      return(HMM)
    }
    HMM$TPM = new_TPM
    HMM$initial_probs = new_initial_probs
    HMM$emis_param$mean = new_mean
    HMM$emis_param$sd = new_sd
  }
}

#' @export
HMM_EM.BernulliHMM <- function(HMM, observations){
  Time = length(observations)
  while(T){
    FB_alg = ForwardBackwardAlgorithm(HMM, observations)
    N_ij = matrix(nrow=HMM$nstates, ncol=HMM$nstates)
    for (i in 1:HMM$nstates){
      for (j in 1:HMM$nstates){
        N_ij[i,j] = with(FB_alg, 
                         sum(Alpha[i,-Time] * Beta[j,-1] * P_mat[j,-1] * HMM$TPM[i,j] / Likelihood_forw))
      }
    }
    N_ia = with(FB_alg,
                rowSums(Alpha * Beta) / Likelihood_forw)
    N_aa = with(FB_alg,
                as.numeric((Alpha * Beta) %*% observations) / Likelihood_forw)
    new_initial_probs = with(FB_alg,
                             (Alpha[,1] * Beta[,1])/Likelihood_forw)
    new_TPM = N_ij / rowSums(N_ij)
    new_p = N_aa / N_ia
    if (sum(abs(new_TPM - HMM$TPM)) <= 0.0001 &&
        sum(abs(new_initial_probs - HMM$initial_probs)) <= 0.0001 &&
        sum(abs(new_lambda - HMM$emis_param$lambda)) <= 0.0001){
      return(HMM)
    }
    HMM$TPM = new_TPM
    HMM$initial_probs = new_initial_probs
    HMM$emis_param$p = new_p
  }
}

