##' psplineRWprior function
##'
##' A function to define Gaussian priors for omega. This function simply stores a vector of means and standard deviations to be passed to the main MCMC function, survspat.
##'
##' @param taumean the prior mean, a vector of length 1 or more. 1 implies a common mean.
##' @param tausd the prior standard deviation, a vector of length 1 or more. 1 implies a common standard deviation.
##' @param basehaz an object inheriting class "basehazardspec", specificlly, this function was used for such objects created by a call to the function PsplineHaz
##' @param order the order of the random walk, default is 2
##' @return an object of class "omegapriorGauss"
##' @seealso \link{survspat}, \link{betapriorGauss}, \link{omegapriorGauss}, \link{etapriorGauss}, \link{indepGaussianprior}, \link{derivindepGaussianprior}
##' @export

psplineRWprior <- function(taumean,tausd,basehaz,order=2){

    warning("!!!UNDER DEVELOPMENT, YOU HAVE BEEN WARNED !!!")

    npar <- distinfo(basehaz)()$npars - 1

    pascal <- choose(order,0:order)
    entries <- suppressWarnings(pascal*c(1,-1)) # alternate signs
    lent <- length(entries)
    d <- t(sapply(0:(npar-lent),function(x){c(rep(0,x),entries,rep(0,(npar-lent)-x))}))

    K <- t(d)%*%d
    #print(K)

    retlist <- list()

    retlist$eval <- function(pars){
        n <- length(pars)
        tau2 <- exp(pars[n])^2
        gamma <- exp(pars[1:(n-1)])
        return(-1/(2*tau2)*gamma%*%K%*%gamma - dnorm(pars[n],taumean,tausd,log=TRUE))
    }

    retlist$deriv <- function(pars){
        n <- length(pars)
        tau <- exp(pars[n])
        gamma <- exp(pars[1:(n-1)])
        deriv1 <- c(gamma*(-1/(2*tau^2)*gamma%*%(K+t(K))),tau*1/(tau^3)*gamma%*%K%*%gamma - 1/tausd^2 * (pars[n]-taumean)) 
        # in the above, gamma = exp(pars[1:(n-1)]) and tau = exp(pars[n]) are the Jacobians
        deriv2 <- as.matrix(Matrix::bdiag(-1/(2*tau^2)*(K+t(K))%*%outer(gamma,gamma),tau*1/(tau^3)*gamma%*%K%*%gamma - 1/tausd^2 * (pars[n]-taumean)) )

        #browser()

        return(list(deriv1=deriv1,deriv2=deriv2))
    }

    return(retlist)

}

##' psplineprior function
##'
##' A function for evaluating the log of an independent Gaussian prior for a given set of parameter values.
##'
##' @param beta parameter beta at which prior is to be evaluated 
##' @param omega parameter omega at which prior is to be evaluated
##' @param eta parameter eta at which prior is to be evaluated
##' @param priors an object of class mcmcPriors, see ?mcmcPriors
##' @return the log of the prior evaluated at the given parameter values
##' @seealso \link{survspat}, \link{betapriorGauss}, \link{omegapriorGauss}, \link{etapriorGauss}, \link{indepGaussianprior}, \link{derivindepGaussianprior}
##' @export

psplineprior <- function(beta=NULL,omega=NULL,eta=NULL,priors){
    
    lp <- 0
    if(!is.null(priors$betaprior)){
        lp <- lp + sum(dnorm(beta,priors$betaprior$mean,priors$betaprior$sd,log=TRUE))
    }

    if(!is.null(priors$omegaprior)){
        #lp <- lp + sum(dnorm(omega,priors$omegaprior$mean,priors$omegaprior$sd,log=TRUE))
        lp <- lp + priors$omegaprior$eval(omega)
    }
    
    if(!is.null(priors$etaprior)){
        lp <- lp + sum(dnorm(eta,priors$etaprior$mean,priors$etaprior$sd,log=TRUE))
    }
    
    return(lp)
}


##' derivpsplineprior function
##'
##' A function for evaluating the first and second derivatives of the log of an independent Gaussian prior
##'
##' @param beta a vector, the parameter beta
##' @param omega a vector, the parameter omega
##' @param eta a vector, the parameter eta 
##' @param priors an object of class 'mcmcPrior', see ?mcmcPrior
##' @return returns the first and second derivatives of the prior
##' @seealso \link{survspat}, \link{betapriorGauss}, \link{omegapriorGauss}, \link{etapriorGauss}, \link{indepGaussianprior}, \link{derivindepGaussianprior}
##' @export

derivpsplineprior <- function(beta=NULL,omega=NULL,eta=NULL,priors){
    if(is.null(beta)){
        deriv1beta <- c()
        deriv2beta <- c()
    }
    else{
        deriv1beta <- (-1/priors$betaprior$sd^2)*(beta-priors$betaprior$mean)
        sdbeta <- priors$betaprior$sd
        if (length(priors$betaprior$sd)<length(beta)){
            sdbeta <- rep(priors$betaprior$sd,length(beta))
        }
        deriv2beta <- -1/sdbeta^2
    }
    if(is.null(omega)){
        deriv1omega <- c()
        deriv2omega <- c()
    }
    else{
        omegastuff <- priors$omegaprior$deriv(omega)
        deriv1omega <- omegastuff$deriv1
        deriv2omega <- omegastuff$deriv2
    }
    if(is.null(eta)){
        deriv1eta <- c()
        deriv2eta <- c()

    }
    else{
        deriv1eta <- (-1/priors$etaprior$sd^2)*(eta-priors$etaprior$mean)
        sdeta <- priors$etaprior$sd
        if (length(priors$etaprior$sd)<length(eta)){
            sdeta <- rep(priors$etaprior$sd,length(eta))
        }
        deriv2eta <- -1/sdeta^2
    }

    deriv1 <- c(deriv1beta,deriv1omega,deriv1eta)
    deriv2 <- diag(c(deriv2beta,deriv2omega,deriv2eta)) # in fact, the 2nd derivative with respect to eta is not necessary, as this is dealt with elsewhere.
    
    return(list(deriv1=deriv1,deriv2=deriv2))
}