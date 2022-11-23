# roxygen2::roxygenise()
# devtools::build(vignettes=TRUE)
# devtools::install(build_vignettes = TRUE)

# ---- ... ----

###Function for checking imports -----------------------------------------------
# Imports <- function() {
#   if (!require(expm)) install.packages('expm')
#   library(expm)
#   if (!require(magrittr)) install.packages('magrittr')
#   library(magrittr)
#   if (!require(dplyr)) install.packages('dplyr')
#   library(dplyr)
#   if (!require(reshape2)) install.packages('reshape2')
#   library(reshape2)
#   if (!require(ggplot2)) install.packages('ggplot2')
#   library(ggplot2)
# }

### HMM Constructor ------------------------------------------------------------
#' @title Hidden Markov Models
#'
#' @description Creates an object of class \code{HMM} containing the parameters
#' of a hidden Markov model with a standard emission distribution.
#' 
#' @usage HMM(distribution, emis_param, initial_probs, TPM)
#' 
#' @param distribution a character of length one specifying the emission distribution,
#' must be one of the following: "Poisson", "Bernoulli", "Normal".
#' @param emis_param the parameters of the emission distribution.
#' For "Poisson": must be a numeric vector or list of length one containing the rates.
#' For "Bernoulli": must be a numeric vector or list of length one containing the probabilities.
#' For "Normal": must be a list with two elements named "mean" and "sd".
#' @param initial_probs a probability vector with the initial distribution of the hidden Markov chain.
#' @param TPM a transition probability matrix for the hidden Markov chain.
#'
#' @details An object of class \code{HMM} also has a subclass, which is either
#' \code{PoissonHMM}, \code{BernoulliHMM}, or \code{NormalHMM} depending on the
#' emission distribution specified in \code{distribution}.
#' 
#' @return \code{HMM} returns an object of class \code{c("PoissonHMM", "HMM")},
#' \code{c("BernoulliHMM", "HMM")}, or \code{c("NormalHMM", "HMM")}.
#' 
#' An object of class \code{HMM} is a list containing the following elements:
#' \describe{
#' \item{\code{nstates}}{a numeric of length one. The size of the hidden state space.}
#' \item{\code{distribution}}{a character of length one. The emission distribution.}
#' \item{\code{emis_param}}{a list containing the parameters of the emission
#' distribution named in accordance with the parameter names used in base R's
#' functions for the probability distributions.}
#' \item{\code{initial_probs}}{a numeric vector of length \code{nstates}. The
#' initial probabilities of the hidden Markov chain.}
#' \item{\code{TPM}}{a square matrix with dimension \code{nstates}. The transition
#' probabilities of the hidden Markov chain.}
#' }
#' 
#' @seealso \code{\link{is.HMM}}, \code{\link{is.PoissonHMM}},
#' \code{\link{is.BernoulliHMM}}, and \code{\link{is.NormalHMM}}
#' for class checking.
#' 
#' @example
#' ## Constructing a two-state Poisson HMM:
#' lambdas <- c(5, 15)
#' delta <- c(0.3, 0.7)
#' P <- rbind(c(0.2, 0.8), c(0.5, 0.5))
#' pois2HMM <- HMM("Poisson", lambdas, delta, P)
#' 
#' @export
HMM <- function(distribution, emis_param, initial_probs, TPM){
  ## Check TPM:
  is.matrix(TPM) || stop("TPM must be a quadratic matrix")
  (dim(TPM)[1] == dim(TPM)[2]) || stop("TPM must be a quadratic matrix")
  all(TPM >= 0) || stop("elements of TPM must be non-negative")
  all(rowSums(TPM) == 1) || stop("rows of TPM must sum to one")
  ## Check initial_probs:
  is.numeric(initial_probs) || stop("initial_probs must be a numeric vector")
  (length(initial_probs) == dim(TPM)[1]) || stop("dimensions of initial_probs and TPM must match")
  all(initial_probs >= 0) || stop("elements of initial_probs must be non-negative")
  (sum(initial_probs) == 1) || stop("initial_prob must be a probability vector")
  
  ## Check distribution:
  is.character(distribution) || stop("distribution must be a character")
  (length(distribution) == 1) || stop("invalid distribution")
  
  l <- list(nstates = length(initial_probs),
            distribution = distribution,
            emis_param = emis_param,
            initial_probs = initial_probs,
            TPM = TPM)
  
  dists <- c("Poisson", "Bernoulli", "Normal")
  sum(distribution == dists) || stop("distribution must be one of following strings: \"Poisson\", \"Bernoulli\", \"Normal\"")
  if (distribution == "Poisson"){
    ## Check parameters:
    if (is.numeric(emis_param)){
      (length(emis_param) == l$nstates) || stop("number of rates in Poisson HMM must correspond to size of hidden state space")
      all(emis_param > 0) || stop("rates for Poisson HMM must be positive")
      l$emis_param <- list(lambda = emis_param)
    } else if (is.list(emis_param)){
      (length(emis_param) == 1) || stop("emis_param does not have the correct number of elements")
      (length(emis_param[[1]]) == l$nstates) || stop("number of rates in Poisson HMM must correspond to size of hidden state space")
      all(unlist(emis_param) > 0) || stop("rates for Poisson HMM must be positive")
      l$emis_param <- list(lambda = emis_param[[1]])
    } else{
      stop("emis_param is not valid for Poisson HMM")
    }
    class(l) <- c("PoissonHMM", "HMM")
  } else if (distribution == "Bernoulli"){
    ## Check parameters:
    if (is.numeric(emis_param)){
      all((emis_param >= 0) && (emis_param <= 1)) || stop("probabilities for Bernoulli HMM must be between 0 and 1")
      (length(emis_param[[1]]) == l$nstates) || stop("number of probabilities in Bernoulli HMM must correspond to size of hidden state space")
      l$emis_param <- list(prob = emis_param)
    } else if (is.list(emis_param)){
      (length(emis_param) == 1) || stop("emis_param does not have the correct number of elements")
      (length(emis_param[[1]]) == l$nstates) || stop("number of probabilities in Bernoulli HMM must correspond to size of hidden state space")
      all((unlist(emis_param) >= 0) && (unlist(emis_param) <= 1)) || stop("probabilities for Bernoulli HMM must be between 0 and 1")
      l$emis_param <- list(prob = emis_param[[1]])
    } else{
      stop("emis_param is not valid for Bernoulli HMM")
    }
    class(l) <- c("BernoulliHMM", "HMM")
  } else if (distribution == "Normal"){
    is.list(emis_param) || stop("emis_param must be a list containing mean and sd for Normal HMM")
    ## Check parameter names:
    if (!(all(names(emis_param) == c("mean", "sd")) || all(names(emis_param) == c("sd", "mean")))){
      stop("emis_param must be a list containing mean and sd for Normal HMM")
    }
    ## Checking parameters:
    all(emis_param$sd > 0) || stop("sd in Normal HMM must be positive")
    ## Change order in list:
    l$emis_param <- list(mean = emis_param$mean, sd = emis_param$sd)
    class(l) <- c("NormalHMM", "HMM")
  }
  return(l)
}


### Function for checking if HMM -----------------------------------------------
#' @title Class Checking for HMM
#' 
#' @description Checks whether an object inherits from class \code{HMM}.
#' 
#' @usage is.HMM(obj)
#' 
#' @seealso \code{\link{HMM}} for object construction, and \code{\link{is.PoissonHMM}},
#' \code{\link{is.BernoulliHMM}}, \code{\link{is.NormalHMM}} for class checking.
#' 
#' @export
is.HMM <- function(obj){
  inherits(obj, "HMM")
}

#' @title Class Checking for PoissonHMM
#' 
#' @description Checks whether an object inherits from class \code{PoissonHMM}.
#' 
#' @usage is.PoissonHMM(obj)
#' 
#' @seealso \code{\link{HMM}} for object construction, and
#' \code{\link{is.BernoulliHMM}}, \code{\link{is.NormalHMM}} for class checking.
#' 
#' @export
is.PoissonHMM <- function(obj){
  inherits(obj, "PoissonHMM")
}

#' @title Class Checking for BernoulliHMM
#' 
#' @description Checks whether an object inherits from class \code{BernoulliHMM}.
#' 
#' @usage is.BernoulliHMM(obj)
#' 
#' @seealso \code{\link{HMM}} for object construction, and \code{\link{is.PoissonHMM}},
#' \code{\link{is.NormalHMM}} for class checking.
#' 
#' @export
is.BernoulliHMM <- function(obj){
  inherits(obj, "BernoulliHMM")
}

#' @title Class Checking for NormalHMM
#' 
#' @description Checks whether an object inherits from class \code{NormalHMM}.
#' 
#' @usage is.NormalHMM(obj)
#' 
#' @seealso \code{\link{HMM}} for object construction, and \code{\link{is.PoissonHMM}},
#' \code{\link{is.BernoulliHMM}}, for class checking.
#' 
#' @export
is.NormalHMM <- function(obj){
  inherits(obj, "NormalHMM")
}



### Print method for HMM--------------------------------------------------------
#' @import dplyr
#' @export
print.HMM <- function(x, ...){
  cat("## Emission distribution:",x$distribution, 
      "\n## Distribution Parameters for emission probability: ",
      unlist(x$emis_param),"\n## Start Probability for the n-states: ",
      unlist(x$initial_probs), "\n## Number of hidden states: ",
      x$nstates, "\n## Class", class(x),"\n")
  cat("## Transition Probability Matrix \n")
  print(x$TPM)
  paste("\n## Probability to go from,",row(x$TPM), "state to",col(x$TPM),
        "state: ",x$TPM) %>% cat()
  cat("\n")
  invisible(x)
}

### More flexible class of HMM taking functions --------------------------------
#' @title Hidden Markov Models
#'
#' @description Creates an object of class \code{xHMM} containing the parameters
#' of a hidden Markov model (HMM) with a user specified emission distribution.
#' 
#' @usage xHMM(distribution, emis_param, emis_paramspace, initial_probs, TPM)
#' 
#' @param distribution a function, the pdf or pmf of the emission distribution.
#' @param emis_param the parameters of the emission distribution.
#' @param emis_paramspace a list containing the parameter spaces
#' of the emission parameters. Elements should be named according to the argument
#' name in \code{distribution} and contain precisely one of the characters:
#' "real", "pos_real", "prob"
#' @param initial_probs a probability vector with the initial distribution of the hidden Markov chain.
#' @param TPM a transition probability matrix for the hidden Markov chain.
#' 
#' @return \code{xHMM} returns an object of class \code{xHMM}.
#' 
#' An object of class \code{xHMM} is a list containing the following elements:
#' \describe{
#' \item{\code{nstates}}{a numeric of length one. The size of the hidden state space.}
#' \item{\code{distribution}}{a function. The emission distribution.}
#' \item{\code{emis_param}}{a list containing the parameters of the emission
#' distribution named in accordance with the parameter names used in the arguments
#' of \code{distribution}.}
#' \item{\code{emis_paramspace}}{a list containing the parameter spaces
#' of the emission parameters. Elements are named according to the argument
#' name in \code{distribution}.}
#' \item{\code{initial_probs}}{a numeric vector of length \code{nstates}. The
#' initial probabilities of the hidden Markov chain.}
#' \item{\code{TPM}}{a square matrix with dimension \code{nstates}. The transition
#' probabilities of the hidden Markov chain.}
#' }
#' 
#' @seealso \code{\link{is.xHMM}} for class checking and \code{\link{HMM}} for
#' Poisson, Bernoulli, or Normal HMMs.
#' 
#' 
#' @export
xHMM <- function(distribution, emis_param, emis_paramspace, initial_probs, TPM){
  ## Check TPM:
  is.matrix(TPM) || stop("TPM must be a quadratic matrix")
  (dim(TPM)[1] == dim(TPM)[2]) || stop("TPM must be a quadratic matrix")
  (all(rowSums(TPM) == 1)) || stop("rows of TPM must sum to one")
  ## Check initial_probs:
  is.numeric(initial_probs) || stop("initial_probs must be a numeric vector")
  (length(initial_probs) == dim(TPM)[1]) || stop("dimensions of initial_probs and TPM must match")
  (sum(initial_probs) == 1) || stop("initial_prob must be a probability vector")
  
  is.function(distribution) || stop("distribution must be a function")
  
  l <- list(nstates = length(initial_probs),
            distribution = distribution,
            emis_param = emis_param,
            emis_paramspace = emis_paramspace,
            initial_probs = initial_probs,
            TPM = TPM)
  class(l) <- c("xHMM")

  return(l)
}

### Class checking ---------------------------
#' @title Class Checking for xHMM
#' 
#' @description Checks whether an object inherits from class \code{xHMM}.
#' 
#' @usage is.xHMM(obj)
#' 
#' @seealso \code{\link{xHMM}} for object construction.
#' 
#' @export
is.xHMM <- function(obj){
  inherits(obj, "xHMM")
}



### Print function for the object xHMM -----------------------------------------
#' @export
print.xHMM <- function(x, ...){
  newstring = substring(deparse(x$distribution)[2], 9)
  cat("## Emission distribution:",gsub("\\,.*","",newstring)) 
  cat("\n## Distribution Parameters for emission probability: ",
      unlist(x$emis_param),"\n## Start Probability for the n-states: ",
      unlist(x$initial_probs), "\n## Number of hidden states: ",x$nstates,
      "\n## The value set of Distribution Parameters : ",
      unlist(x$emis_paramspace), "\n## Class", class(x),"\n")
  cat("\n## Transition Probability Matrix \n")
  print(x$TPM)
  paste("\n## Probability to go from,",row(x$TPM), "state to",col(x$TPM),
        "state: ",x$TPM) %>% cat()
  cat("\n")
  invisible(x)
}

