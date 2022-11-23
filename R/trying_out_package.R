# This document is only for "testing" purposes
library(magrittr)
# Testing specified Poisson HMM ----
pois_example = HMM(distribution = "Poisson",
                   emis_param = list(lambda=c(10, 20)),
                   initial_probs = rep(1/2, 2),
                   TPM = matrix(1/2, 2, 2))
pois_example

earthquakes = read.table("http://www.hmms-for-time-series.de/second/data/earthquakes.txt")[,2]

pois_fit = HMM_EM(HMM = pois_example, observations = earthquakes)

pois_fit

ForwardBackwardAlgorithm(pois_fit, earthquakes)$Likelihood_back

# Testing more flexible HMM and general nlm method ----
pois_example_2 = xHMM(distribution = dnorm,
                      emis_param = list(mean=c(10, 20),sd=c(2,6)),
                      emis_paramspace = list(mean="real",sd="pos_real"),
                      initial_probs = rep(1/2, 2),
                      TPM = matrix(1/2, 2, 2))

class(pois_example_2)

ForwardBackwardAlgorithm(pois_example_2,earthquakes)

pois_example_2
pois_fit_nlm = nlm_hmm(HMM = pois_example_2,
                       observations = earthquakes)

pois_fit_nlm
# Comparing the results of the two ----
pois_fit_nlm$TPM
pois_fit$TPM
pois_fit_nlm$initial_probs
pois_fit$initial_probs
pois_fit$emis_param
pois_fit_nlm$emis_param

rpoissonHMM = function(n,HMM){
  X = c(1:n)
  State = c(1:n)
  State[1] = sample(c(1,2),1,prob = HMM$initial_probs)
  X[1] = rpois(1,lambda = HMM$emis_param$lambda[State[1]])
  if (n==1){return(X)}
  for (i in c(2:n)){
    State[i] =sample(c(1,2),1,prob = HMM$TPM[,State[i-1]])
    X[i] = rpois(1,lambda = HMM$emis_param$lambda[State[i]])
  }
  return(X)
}
l = set.seed(1)
data = rpoissonHMM(100,pois_example)
alpha = ForwardBackwardAlgorithm(pois_example, data)

#alpha = alpha$Likelihood_forw
alpha$LogLikelihood_forw
alpha$LogLikelihood_back
alpha1 = ForwardBackwardAlgorithm(pois_example, data)$Alpha
log(sum(alpha1[,dim(alpha1)[2]]))

Local_decoding(pois_example,earthquakes)
Local_decoding(pois_example_2,earthquakes)

Forecasting(pois_example,earthquakes)
Forecasting(pois_example_2,earthquakes)

State_prediction(pois_example,earthquakes,step=100)
State_prediction(pois_example_2,earthquakes, step = 1000)
# Testing nstates_HMM function to for easy fitting of data to n states:

test = Nstate_HMM(distribution = "Poisson",
           data = earthquakes,
           nstates = 3)
test$TPM

# Find AICBIC for states up to n

Nstates_AICBIC("Poisson", earthquakes, 7, "base")
Nstates_AICBIC("Poisson", earthquakes, 7, "ggplot")

