#' @import ggplot2 reshape2 dplyr


#' @title AIC and BIC of Specified HMM
#' @description calculates the AIC and BIC of class HMM or xHMM.
#' 
#' 
#' @usage AICBIC(HMM, observations)
#' 
#' @param HMM An object of class \code{HMM} or \code{xHMM}, a hidden Markov model.
#' @param observations A numeric vector containing the observed data.
#'
#' @return A list containing \code{AIC} and \code{BIC}.
#' @export
AICBIC <- function(HMM,observations){
  log_L = ForwardBackwardAlgorithm(HMM, observations)$LogLikelihood_forw
  n_obs = length(observations)
  nstates = length(HMM$initial_probs)
  p = nstates*(nstates-1) + nstates*length(HMM$emis_param)
  return(list(AIC = -2*log_L+2*p,
              BIC = -2*log_L+p*log(n_obs)))
}

#' @title Find Stationary Distribution of HMM
#' 
#' @usage Stationary(HMM)
#'
#' @param HMM the \code{HMM} or \code{xHMM} object to determine the stationary
#' distribution of.
#'
#' @return a \code{HMM} or \code{xHMM} object identical to the input except the
#' \code{initial_probs} parameter has been changed to the stationary distribution
#' of the (hidden) Markov chain.
#' @export
#'
Stationary <- function(HMM) {
  transition = HMM$TPM
  m = dim(transition)[1]
  res = rep(1, m) %*% solve(diag(m) - transition + matrix(1, m, m))
  HMM$initial_probs = res
  return(HMM)
}


#' @title Calculating AIC and BIC for Different State Space Sizes
#' 
#' @usage Nstates_AICBIC(distribution = "Normal", data, nstates = 3, doplot = "ggplot")
#'
#' @param distribution a character of length one specifying the emission distribution,
#' must be one of the following: "Poisson", "Bernoulli", "Normal".
#' @param data a numeric vector containing the observed data.
#' @param nstates the number of states up to which the function should fit using the function \code{\link{nstates_HMM}}.
#' i.e. for \code{nstates=5} the data is fit and the AIC and BIC calculated for states 2 up to 5.
#' @param doplot can be specified as either "ggplot" or "base" according to which plot is wanted.
#' otherwise no plot will be made.
#'
#' @return Returns a \code{data.frame} containing AIC and BIC for HMM with the
#' specified \code{distribution} and stationary initial distribution for the
#' hidden Markov chain, when fitted to the data for 2 states up to \code{nstates}
#' states. Furthermore a plot is created in accordance with what has been specified
#' by \code{doplot}.
#' @export
#'
Nstates_AICBIC <- function(distribution="Normal", data, nstates=3, doplot="ggplot"){
  AICBIC = data.frame(matrix(nrow=0, ncol=2))
  colnames(AICBIC) = c("AIC", "BIC")
  for (i in 2:nstates){
    AICBIC = rbind(AICBIC, AICBIC(Stationary(Nstate_HMM(distribution, data, i)), data))
  }
  AICBIC = cbind("NoStates"=2:nstates, AICBIC)
  if (doplot=="base"){
    plot(rep(AICBIC$"NoStates", 2), c(AICBIC$AIC, AICBIC$BIC),
         col=c(rep('red', nstates-1), rep('blue', nstates-1)),
         xlab="Number of States", ylab="Information Criterion Value",
         main=c("AIC and BIC", paste(distribution, "Distribution; 2 to", nstates, "states")))
    lines(2:nstates, AICBIC$AIC, col="red")
    lines(2:nstates, AICBIC$BIC, col="blue")
    legend("topleft", inset=.02, c("AIC","BIC"), fill=c("red", "blue"), horiz=TRUE)
  }
  if (doplot=="ggplot"){
    print(ggplot2::ggplot(reshape2::melt(AICBIC, id="NoStates"), ggplot2::aes(x=NoStates, y=value, color=variable)) 
          + ggplot2::geom_point() + ggplot2::geom_line() + ggplot2::ggtitle(paste("AIC and BIC for", distribution, "distribution; 2 to", nstates, "states.")))
  }
  return(AICBIC)
}
