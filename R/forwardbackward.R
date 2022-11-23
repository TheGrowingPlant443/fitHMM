forwardbackwardalgorithm <- function(HMM, observations){
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