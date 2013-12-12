##' mcmcPriors function
##'
##' A function to define priors for mcmc
##'
##' @param betaprior prior for beta, the covariate effects
##' @param omegaprior prior for omega, the parameters of the baseline hazard
##' @param etaprior prior for eta, the parameters of the latent field
##' @param call function to evaluate the log-prior e.g. logindepGaussianprior
##' @param derivative function to evaluate the first and second derivatives of the prior 
##' @return an onject of class mcmcPriors
##' @seealso \link{survspat},
##' @export

mcmcPriors <- function(betaprior=NULL,omegaprior=NULL,etaprior=NULL,call=NULL,derivative=NULL){
    retlist <- list()
    retlist$betaprior <- betaprior
    retlist$omegaprior <- omegaprior
    retlist$etaprior <- etaprior
    retlist$call <- call
    retlist$derivative <- derivative
    class(retlist) <- "mcmcPriors"
    return(retlist)
}



##' betapriorGauss function
##'
##' A function to define Gaussian priors for beta.
##'
##' @param mean the prior mean, a vector of length 1 or more. 1 implies a common mean.
##' @param sd the prior standard deviation, a vector of length 1 or more. 1 implies a common standard deviation.
##' @return an object of class "betapriorGauss"
##' @export

betapriorGauss <- function(mean,sd){
    retlist <- list()
    retlist$mean <- mean
    retlist$sd <- sd
    class(retlist) <- "betapriorGauss"
    return(retlist) 
}



##' omegapriorGauss function
##'
##' A function to define Gaussian priors for omega.
##'
##' @param mean the prior mean, a vector of length 1 or more. 1 implies a common mean.
##' @param sd the prior standard deviation, a vector of length 1 or more. 1 implies a common standard deviation.
##' @return an object of class "omegapriorGauss"
##' @export

omegapriorGauss <- function(mean,sd){
    retlist <- list()
    retlist$mean <- mean
    retlist$sd <- sd
    class(retlist) <- "omegapriorGauss"
    return(retlist) 
}



##' etapriorGauss function
##'
##' A function to define Gaussian priors for eta.
##'
##' @param mean the prior mean, a vector of length 1 or more. 1 implies a common mean.
##' @param sd the prior standard deviation, a vector of length 1 or more. 1 implies a common standard deviation.
##' @return an object of class "etapriorGauss"
##' @export

etapriorGauss <- function(mean,sd){
    retlist <- list()
    retlist$mean <- mean
    retlist$sd <- sd
    class(retlist) <- "etapriorGauss"
    return(retlist) 
}