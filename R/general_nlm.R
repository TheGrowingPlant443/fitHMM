# ---- Transform Parameters ----

transform_param <- function(parameters, transformation){
  transformed = parameters
  parameter_names = names(parameters)
  for (n in parameter_names){
    if (transformation[n] == "real"){
      # Parameters  should remain the same
    }
    else if (transformation[n] == "pos_real"){
      # Transform using log:
      transformed[[n]] = log(parameters[[n]])
    }
    else if (transformation[n] == "prob"){
      # Transform dividing by first entry, then taking log
      transformed[[n]] = log(-(log(parameters[[n]])-1e-100)) # makes sure no entries are -inf
    }
    else{
      stop("Problem in specification of parameters or their transformation")
    }
  }
  return(transformed)
}

detransform_param <- function(parameters, transformation){
  detransformed = transformation
  parameter_names = names(transformation)
  for (n in parameter_names){
    if (transformation[n] == "real"){
      # Parameters  should remain the same
      detransformed[[n]] = parameters[[n]]
    }
    else if (transformation[n] == "pos_real"){
      # Transform using log:
      detransformed[[n]] = exp(parameters[[n]])
    }
    else if (transformation[n] == "prob"){
      # Transform dividing by first entry, then taking log
      detransformed[[n]] = exp(-exp(parameters[[n]]))
    }
    else{
      stop("Problem in specification of parameters or their transformation")
    }
  }
  return(detransformed)
}

# ---- Vectorfy and Devectorfy ----

vectorfy <- function(delta, TPM, param, trans){
  n_param = length(param[[1]])
  m = length(delta)
  my_vector = log(delta[-m]/delta[m])
  my_vector = append(my_vector, log(TPM[!diag(m)]/diag(TPM)))
  transformed_params = transform_param(param, trans)
  for (param in transformed_params){
    my_vector = append(my_vector, param)
  }
  return(my_vector)
}

devectorfy <- function(my_vector, n_param, m, Time, trans){
  param_names = names(trans)
  delta = c(exp(my_vector[1:(m-1)]), 1)
  delta = delta/sum(delta)
  TPM = diag(m)
  TPM[!TPM] = exp(my_vector[m:(m+(m-1)*m-1)])
  TPM = TPM/apply(TPM, 1, sum)
  parms = list()
  for (i in 1:m){
    parms[[param_names[i]]] = my_vector[(((m+(m-1)*m-1)+1)+(i-1)*n_param):((m+(m-1)*m-1)+i*n_param)]
  }
  parms = detransform_param(parms, trans)
  return(list(delta = delta,
              TPM = TPM,
              parameters = parms))
}

# ---- Forward Algorithm - log-sum-exp procedure ----

ForwardAlgorithm_nlm <- function(dist, param, obs, TPM, delta){
  S = c()
  m = length(delta)
  Time = length(obs)
  P_mat = matrix(nrow=m, ncol=Time)
  for (i in 1:Time){
    P_mat[,i] = do.call(dist, append(list(obs[i]), param))
  }
  Alpha = matrix(nrow=m, ncol=Time)
  Alpha[,1] = matrix(delta * P_mat[,1])
  S[1] = sum(Alpha[,1])
  for (t in 2:Time){
    Alpha[,t] = (Alpha[,t-1] %*% TPM * P_mat[,t])/S[t-1]
    S[t] = sum(Alpha[,t])
  }
  return(list(Alpha_Matrix = Alpha,
              LogLikelihood = sum(log(S))))
}

#' @title Numeric Parameter Estimation of HMM
#' 
#' @description Uses \code{nlm} to minimize the negative log likelihood of an HMM,
#' which corresponds to find the parameters maximizing the likelihood.
#' 
#' @usage nlm_hmm(HMM, observations)
#'
#' @param HMM object of class \code{xHMM} containing the initial parameters.
#' @param observations numeric vector containing the observed sequence.
#'
#' @return Returns a new object of class \code{xHMM} that is fitted to the
#' observations according to non-linear minimization approach.
#' @export
#'
nlm_hmm <- function(HMM, observations){
  # some new references to the information in the HMM object is made for ease of notation:
  dist = HMM$distribution
  param = HMM$emis_param
  param_trans = HMM$emis_paramspace
  TPM = HMM$TPM
  delta = HMM$initial_probs
  obs = observations
  # The likelihood function written in a way such that nlm() can find the desired parameters:
  likelihood <- function(parameters, observations){
    # "parameters" is a numeric vector containing all the parameters:
    # In order to use the forward algorithm these need to be split into the correct form, which is done by the function devectorfy:
    devectorfied = devectorfy(parameters, n_param, m, Time, param_trans)
    # devectorfy returns a list containing the desired information:
    TPM_t = devectorfied$TPM
    delta_t = devectorfied$delta
    parameters_t = devectorfied$parameters
    LogxHMM = xHMM(dist,param, param_trans, delta, TPM)
    LogxHMM$emis_param = parameters_t
    LogxHMM$TPM = TPM_t
    LogxHMM$initial_probs = delta_t
    # the negative log likelihood is returned using the given parameters and observations:
    return(-ForwardAlgorithm_nlm(dist, parameters_t, observations, TPM_t, delta_t)$LogLikelihood)
    #return(-ForwardBackwardAlgorithm(LogxHMM,observations)$LogLikelihood_forw)
  }
  # Find the number of parameters for each state
  n_param = length(param[[1]])
  # finding the number of states m and the number of observations Time:
  m = length(delta)
  Time = length(obs)
  # In order to use nlm to optimize all the parameters delta, TPM, and the state parameters are combined into a single numeric vector:
  my_vector = vectorfy(delta, TPM, param, param_trans)
  # The parameters are optimized using nlm:
  estimates = nlm(likelihood, my_vector, obs)$estimate
  # Converting the numeric vector containing the parameters into fitted: TPM, emis_param, and initial_probs:
  real_estimates = devectorfy(estimates, n_param, m, Time, param_trans)
  # The original HMM is then updated with the fit, and the HMM is returned:
  # (remark that this makes a copy of HMM and does not change it in the global environment)
  HMM$emis_param = real_estimates$param
  HMM$TPM = real_estimates$TPM
  HMM$initial_probs = real_estimates$delta
  return(HMM)
}

