##' exp_ltar function
##'
##' A function to evaluate the likelihood for the parametric proportional hazards model with baseline hazard derived from the exponential model.
##'
##' @param t vector of times 
##' @param XbetaplusY the matrix product X times beta where  is the design mastrix and beta are the covariate effects
##' @param expXbetaplusY the exponential of the matrix product X times beta where  is the design mastrix and beta are the covariate effects 
##' @param theta rate parameter for exponential model
##' @return the log-likelihood
##' @export

exp_ltar <- function(t,XbetaplusY,expXbetaplusY,theta){
    n <- length(t)
    return(XbetaplusY + log(theta)-expXbetaplusY*theta*t)
}


##' weibull_ltar function
##'
##' A function to evaluate the likelihood for the parametric proportional hazards model with baseline hazard derived from the weibull model.
##'
##' @param t vector of times 
##' @param XbetaplusY the matrix product X times beta where  is the design mastrix and beta are the covariate effects
##' @param expXbetaplusY the exponential of the matrix product X times beta where  is the design mastrix and beta are the covariate effects 
##' @param theta vector of length 2: shape and scale parameters for the weibull model
##' @return the log-likelihood
##' @export

weibull_ltar <- function(t,XbetaplusY,expXbetaplusY,theta){
    n <- length(t)
    alpha <- theta[1]
    lambda <- (1/theta[2])^theta[1]
    return(XbetaplusY + log(lambda) + log(alpha) + t^(alpha-1) - expXbetaplusY*lambda*t^alpha)
}

##' logindepnormalprior function
##'
##' A function to evaluate the log prior for independent normals
##'
##' @param beta parameter beta at which prior is to be evaluated 
##' @param omega parameter omega at which prior is to be evaluated
##' @param betapriormean prior mean for beta 
##' @param betapriorsd prior standard deviation for beta
##' @param omegapriormean prior mean fpr omega 
##' @param omegapriorsd prior standard deviation for omega 
##' @return the log prior
##' @export

logindepnormalprior <- function(beta,omega,betapriormean,betapriorsd,omegapriormean,omegapriorsd){
    return(sum(dnorm(beta,betapriormean,betapriorsd,log=TRUE))+sum(dnorm(omega,omegapriormean,omegapriorsd,log=TRUE)))
}



##' logindepGaussianprior function
##'
##' A function to evaluate the log prior for independent normals
##'
##' @param beta parameter beta at which prior is to be evaluated 
##' @param omega parameter omega at which prior is to be evaluated
##' @param eta parameter eta at which prior is to be evaluated
##' @param priors an object of class mcmcPriors, see ?mcmcPriors
##' @return the log prior
##' @export

logindepGaussianprior <- function(beta=NULL,omega=NULL,eta=NULL,priors){
    
    lp <- 0
    if(!is.null(priors$betaprior)){
        lp <- lp + sum(dnorm(beta,priors$betaprior$mean,priors$betaprior$sd,log=TRUE))
    }

    if(!is.null(priors$omegaprior)){
        lp <- lp + sum(dnorm(omega,priors$omegaprior$mean,priors$omegaprior$sd,log=TRUE))
    }
    
    if(!is.null(priors$etaprior)){
        lp <- lp + sum(dnorm(eta,priors$etaprior$mean,priors$etaprior$sd,log=TRUE))
    }
    
    return(lp)
}



##' logposterior_exponential_nospat function
##'
##' A function to evaluate the log-posterior
##'
##' @param tm vector of observed times
##' @param delta censoring indicator
##' @param X design matrix
##' @param beta beta vector at which to evaluate the posterior
##' @param omega omega vector at which to evaluate the posterior
##' @param betapriormean prior mean for beta 
##' @param betapriorsd prior standard deviation for beta
##' @param omegapriormean prior mean fpr omega 
##' @param omegapriorsd prior standard deviation for omega 
##' @return the log posterior
##' @export

logposterior_exponential_nospat <- function(tm,delta,X,beta,omega,betapriormean,betapriorsd,omegapriormean,omegapriorsd){
    Xbeta <- X%*%beta
    theta <- exp(omega)
    return(sum(delta*(Xbeta+omega)-exp(Xbeta)*theta*tm) + logindepnormalprior(beta=beta,omega=omega,betapriormean=betapriormean,betapriorsd=betapriorsd,omegapriormean=omegapriormean,omegapriorsd=omegapriorsd))
}









