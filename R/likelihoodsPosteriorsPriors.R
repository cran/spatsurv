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




##' logposterior.exp function
##'
##' A function to evaluate the log-posterior
##'
##' @param tm vector of observed times
##' @param delta censoring indicator
##' @param X design matrix
##' @param beta beta vector at which to evaluate the posterior
##' @param omega omegaa vector at which to evaluate the posterior
##' @param eta eta vector at which to evaluate the posterior
##' @param gamma the transformed latent Gaussian field at which to evaluate the posterior
##' @param priors the priors, an object of class mcmcPriors
##' @param covmodel an object of class covmodel, see ?covmodel
##' @param u distance matrix
##' @param control list of control parameters, see ?inference.control
##' @return the log posterior
##' @return ...
##' @export

logposterior.exp <- function(tm,delta,X,beta,omega,eta,gamma,priors,covmodel,u,control){

    n <- nrow(X)
    Xbeta <- X%*%beta
    sigma <- matrix(getcov(u=u,sigma=exp(eta[1]),phi=exp(eta[2]),model=covmodel$model,pars=covmodel$pars),n,n)
    cholsigma <- t(chol(sigma))
    priorcontrib <- -(1/2)*sum(gamma^2) + do.call(priors$call,args=list(beta=beta,omega=omega,eta=eta,priors=priors))
    Y <- -eta[1]^2/2 + cholsigma%*%gamma
    stuff <- Xbeta + Y + omega
    expstuff <- exp(stuff)

    logpost <- sum(delta*(stuff)-expstuff*tm) + priorcontrib
    
    grad <- rep(0,length(beta)+length(omega)+length(eta)+length(gamma))

    deriv <- do.call(priors$derivative,args=list(beta=beta,omega=omega,eta=eta,priors=priors))$deriv1 # first derivaties of priors
    deriv <- c(deriv,-gamma) # tag on gamma
    deriv[(length(beta)+length(omega)+1):(length(beta)+length(omega)+length(eta))] <- 0 # random walk for eta ...
    
    stuff2 <- delta-expstuff*tm
    grad <- rep(0,length(deriv))    
    for(i in 1:length(beta)){
        grad[i] <- grad[i] + sum(X[,i]*stuff2)
    }
    for (i in 1:length(omega)){
        grad[length(beta)+i] <- grad[length(beta)+i] + sum(stuff2)
    }
    for(i in 1:length(gamma)){
        grad[length(beta)+length(omega)+length(eta)+i] <- grad[length(beta)+length(omega)+length(eta)+i] + sum(cholsigma[,i]*stuff2)
    }   
    
    grad <- grad + deriv

    return(list(logpost=logpost,grad=grad,Y=Y))
}


##' logposterior.exp.gridded function
##'
##' A function to evaluate the log-posterior for the gridded model
##'
##' @param tm vector of observed times
##' @param delta censoring indicator
##' @param X design matrix
##' @param beta beta vector at which to evaluate the posterior
##' @param omega omegaa vector at which to evaluate the posterior
##' @param eta eta vector at which to evaluate the posterior
##' @param gamma the transformed latent Gaussian field at which to evaluate the posterior
##' @param priors the priors, an object of class mcmcPriors
##' @param covmodel an object of class covmodel, see ?covmodel
##' @param u distance matrix
##' @param control list of control parameters, see ?inference.control
##' @return the log posterior
##' @return ...
##' @export

logposterior.exp.gridded <- function(tm,delta,X,beta,omega,eta,gamma,priors,covmodel,u,control){

    n <- nrow(X)
    Xbeta <- X%*%beta

    covbase <- matrix(getcov(u=u,sigma=exp(eta[1]),phi=exp(eta[2]),model=covmodel$model,pars=covmodel$pars),control$Mext,control$Next)
    invrootQeigs <- sqrt(Re(fft(covbase)))
    
    Ygrid <- YfromGamma(gamma,invrootQeigs=invrootQeigs,mu=-(exp(eta[1]))^2/2)   
                      
    priorcontrib <- -(1/2)*sum(gamma^2) + do.call(priors$call,args=list(beta=beta,omega=omega,eta=eta,priors=priors))
    stuff <- Xbeta + Ygrid[control$idx] + omega
    expstuff <- exp(stuff)
    logpost <- sum(delta*(stuff)-expstuff*tm) + priorcontrib 
   
    grad <- rep(0,length(beta)+length(omega)+length(eta)+length(gamma))

    deriv <- do.call(priors$derivative,args=list(beta=beta,omega=omega,eta=eta,priors=priors))$deriv1 # first derivaties of priors
    deriv <- c(deriv,-gamma) # tag on gamma
    deriv[(length(beta)+length(omega)+1):(length(beta)+length(omega)+length(eta))] <- 0 # random walk for eta ...
    
    stuff2 <- delta-expstuff*tm
    grad <- rep(0,length(deriv))    
    for(i in 1:length(beta)){
        grad[i] <- grad[i] + sum(X[,i]*stuff2)
    }
    for (i in 1:length(omega)){
        grad[length(beta)+i] <- grad[length(beta)+i] + sum(stuff2)
    }
    
    lenbeta <- length(beta)
    lenomega <- length(omega)
    leneta <- 2
    lenY <- length(Ygrid)
    npars <- lenbeta + lenomega + leneta + lenY
        
    bitsnbobs <- matrix(0,control$Mext,control$Next)
    bitsnbobs[control$uqidx] <- bitsnbobs[control$uqidx] + sapply(control$uqidx,function(i){sum(stuff2[control$idx==i])})
    grad[(lenbeta+lenomega+leneta+1):npars] <- grad[(lenbeta+lenomega+leneta+1):npars] + Re((1/(control$Mext*control$Next))*fft(invrootQeigs*fft(bitsnbobs,inverse=TRUE)))
      
    grad <- grad + deriv

    return(list(logpost=logpost,grad=grad,Y=Ygrid))
}



##' logposterior.weibull function
##'
##' A function to evaluate the log-posterior
##'
##' @param tm vector of observed times
##' @param delta censoring indicator
##' @param X design matrix
##' @param beta beta vector at which to evaluate the posterior
##' @param omega omegaa vector at which to evaluate the posterior
##' @param eta eta vector at which to evaluate the posterior
##' @param gamma the transformed latent Gaussian field at which to evaluate the posterior
##' @param priors the priors, an object of class mcmcPriors
##' @param covmodel an object of class covmodel, see ?covmodel
##' @param u distance matrix
##' @param control list of control parameters, see ?inference.control
##' @return the log posterior
##' @return ...
##' @export

logposterior.weibull <- function(tm,delta,X,beta,omega,eta,gamma,priors,covmodel,u,control){

    alpha <- exp(omega[1])
    lambda <- exp(omega[2])

    n <- nrow(X)
    Xbeta <- X%*%beta
    sigma <- matrix(getcov(u=u,sigma=exp(eta[1]),phi=exp(eta[2]),model=covmodel$model,pars=covmodel$pars),n,n)
    cholsigma <- t(chol(sigma))
    priorcontrib <- -(1/2)*sum(gamma^2) + do.call(priors$call,args=list(beta=beta,omega=omega,eta=eta,priors=priors))
    Y <- -eta[1]^2/2 + cholsigma%*%gamma
    stuff <- Xbeta + Y
    expstuff <- exp(stuff)

    logpost <- sum(delta*(stuff + log(lambda) + log(alpha) + (alpha-1)*log(tm))-expstuff*lambda*tm^alpha) + priorcontrib
    
    grad <- rep(0,length(beta)+length(omega)+length(eta)+length(gamma))

    deriv <- do.call(priors$derivative,args=list(beta=beta,omega=omega,eta=eta,priors=priors))$deriv1 # first derivaties of priors
    deriv <- c(deriv,-gamma) # tag on gamma
    deriv[(length(beta)+length(omega)+1):(length(beta)+length(omega)+length(eta))] <- 0 # random walk for eta ...
    
    stuff2 <- delta-expstuff*lambda*tm^alpha
    grad <- rep(0,length(deriv))    
    for(i in 1:length(beta)){
        grad[i] <- grad[i] + sum(X[,i]*stuff2)
    }
    for (i in 1:length(omega)){
        if(i ==1){
            grad[length(beta)+i] <- grad[length(beta)+i] + sum(delta*(1/alpha+log(tm))-expstuff*lambda*log(tm)*tm^alpha)*alpha # alpha term at end present by chain rule: since alpha = exp(omega[1]) we have dalpha/domega[1] = exp(omega[1]) = alpha  
        }
        else{
            grad[length(beta)+i] <- grad[length(beta)+i] + sum(delta/lambda-expstuff*tm^alpha)*lambda # lambda term at end present by chain rule: since lambda = exp(omega[2]) we have dlambda/domega[2] = exp(omega[2]) = lambda
        }
    }
    for(i in 1:length(gamma)){
        grad[length(beta)+length(omega)+length(eta)+i] <- grad[length(beta)+length(omega)+length(eta)+i] + sum(cholsigma[,i]*stuff2)
    }   
    
    #browser()
    
    grad <- grad + deriv

    return(list(logpost=logpost,grad=grad,Y=Y))
}


##' logposterior.weibull.gridded function
##'
##' A function to evaluate the log-posterior for gridded data
##'
##' @param tm vector of observed times
##' @param delta censoring indicator
##' @param X design matrix
##' @param beta beta vector at which to evaluate the posterior
##' @param omega omegaa vector at which to evaluate the posterior
##' @param eta eta vector at which to evaluate the posterior
##' @param gamma the transformed latent Gaussian field at which to evaluate the posterior
##' @param priors the priors, an object of class mcmcPriors
##' @param covmodel an object of class covmodel, see ?covmodel
##' @param u distance matrix
##' @param control list of control parameters, see ?inference.control
##' @return the log posterior
##' @return ...
##' @export

logposterior.weibull.gridded <- function(tm,delta,X,beta,omega,eta,gamma,priors,covmodel,u,control){

    alpha <- exp(omega[1])
    lambda <- exp(omega[2])

    n <- nrow(X)
    Xbeta <- X%*%beta

    covbase <- matrix(getcov(u=u,sigma=exp(eta[1]),phi=exp(eta[2]),model=covmodel$model,pars=covmodel$pars),control$Mext,control$Next)
    invrootQeigs <- sqrt(Re(fft(covbase)))
    
    Ygrid <- YfromGamma(gamma,invrootQeigs=invrootQeigs,mu=-(exp(eta[1]))^2/2)   
                      
    priorcontrib <- -(1/2)*sum(gamma^2) + do.call(priors$call,args=list(beta=beta,omega=omega,eta=eta,priors=priors))
    stuff <- Xbeta + Ygrid[control$idx]
    expstuff <- exp(stuff)

    logpost <- sum(delta*(stuff + log(lambda) + log(alpha) + (alpha-1)*log(tm))-expstuff*lambda*tm^alpha) + priorcontrib
    
    grad <- rep(0,length(beta)+length(omega)+length(eta)+length(gamma))

    deriv <- do.call(priors$derivative,args=list(beta=beta,omega=omega,eta=eta,priors=priors))$deriv1 # first derivaties of priors
    deriv <- c(deriv,-gamma) # tag on gamma
    deriv[(length(beta)+length(omega)+1):(length(beta)+length(omega)+length(eta))] <- 0 # random walk for eta ...
    
    stuff2 <- delta-expstuff*lambda*tm^alpha
    grad <- rep(0,length(deriv))    
    for(i in 1:length(beta)){
        grad[i] <- grad[i] + sum(X[,i]*stuff2)
    }
    for (i in 1:length(omega)){
        if(i ==1){
            grad[length(beta)+i] <- grad[length(beta)+i] + sum(delta*(1/alpha+log(tm))-expstuff*lambda*log(tm)*tm^alpha)*alpha # alpha term at end present by chain rule: since alpha = exp(omega[1]) we have dalpha/domega[1] = exp(omega[1]) = alpha  
        }
        else{
            grad[length(beta)+i] <- grad[length(beta)+i] + sum(delta/lambda-expstuff*tm^alpha)*lambda # lambda term at end present by chain rule: since lambda = exp(omega[2]) we have dlambda/domega[2] = exp(omega[2]) = lambda
        }
    }
    
    lenbeta <- length(beta)
    lenomega <- length(omega)
    leneta <- 2
    lenY <- length(Ygrid)
    npars <- lenbeta + lenomega + leneta + lenY
        
    bitsnbobs <- matrix(0,control$Mext,control$Next)
    bitsnbobs[control$uqidx] <- bitsnbobs[control$uqidx] + sapply(control$uqidx,function(i){sum(stuff2[control$idx==i])})
    grad[(lenbeta+lenomega+leneta+1):npars] <- grad[(lenbeta+lenomega+leneta+1):npars] + Re((1/(control$Mext*control$Next))*fft(invrootQeigs*fft(bitsnbobs,inverse=TRUE)))
    
    grad <- grad + deriv

    return(list(logpost=logpost,grad=grad,Y=Ygrid))
}




##' logposterior.gamma function
##'
##' 
##'
##' @param tm an gamma object
##' @param ... additional arguments
##' @return ...
##' @export

logposterior.gamma <- function(tm,...){
    stop("not implemented as yet ...")
}
