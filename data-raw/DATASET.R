## code to prepare `DATASET` dataset goes here

usethis::use_data(DATASET, overwrite = TRUE)

## Simulate from a two-state Poisson HMM:
HMM <- list(nstates = 2,
            distribution = "Poisson",
            emis_param = list(lambda = c(10, 20)),
            initial_probs = c(0.5, 0.5),
            TPM = rbind(c(0.6, 0.4), c(0.2, 0.8)))

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

set.seed(1)
Poisson2State <- rpoissonHMM(100, HMM)

usethis::use_data(Poisson2State)

rnbinom_2hmm <- function(n, delta, TPM, sizes, probs){
  observed = c()
  state = rbinom(1, 1, delta[1])
  for (i in 1:n){
    if (state){
      observed[i] = rnbinom(1, sizes[1], probs[1])
      state = rbinom(1, 1, TPM[1,1])
    }
    else{
      observed[i] = rnbinom(1, sizes[2], probs[2])
      state = rbinom(1, 1, TPM[2,1])
    }
  }
  return(observed)
}

Nbinom2State <- rnbinom_2hmm(100, 
                             delta=c(0.5, 0.5), 
                             TPM=rbind(c(0.6, 0.4), c(0.2, 0.8)),
                             sizes=c(10, 100),
                             probs=c(0.25, 0.30))

usethis::use_data(Nbinom2State)
