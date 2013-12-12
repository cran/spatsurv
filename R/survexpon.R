##' survexpon function
##'
##' A function to run a Bayesian analysis on right censored survial data assuming a proportional hazards model with
##' baseline hazard derived from the exponential model.
##'
##' @param formula see ?flexsurvreg 
##' @param data X see ?flexsurvreg
##' @param mcmc.control mcmc control parameters, see ?mcmcpars
##' @param betapriormean prior mean for beta, default is 0
##' @param betapriorsd prior standard deviation for beta
##' @param omegapriormean prior mean for omega, default is 0
##' @param omegapriorsd prior standard deviation for omega 
##' @return the mcmc output
##' @export


survexpon <- function(formula,data,mcmc.control,betapriormean=0,betapriorsd,omegapriormean=0,omegapriorsd){
                        
    ##########
    # This chunk of code borrowed from flexsurvreg    
    ##########                    
    call <- match.call()
    indx <- match(c("formula", "data"), names(call), nomatch = 0)
    if (indx[1] == 0){ 
        stop("A \"formula\" argument is required")
    }
    temp <- call[c(1, indx)]
    temp[[1]] <- as.name("model.frame")
    m <- eval(temp, parent.frame())
    Y <- model.extract(m, "response")
    if (!inherits(Y, "Surv")){ 
        stop("Response must be a survival object")
    }
    if (!(attr(Y, "type") %in% c("right", "counting"))){ 
        stop("Survival object type \"", attr(Y, "type"), "\"", 
            " not supported")
    }
    if (attr(Y, "type") == "counting"){ 
        Y <- cbind(Y, time = Y[, "stop"] - Y[, "start"])
    }
    else{ 
        Y <- cbind(Y, start = 0, stop = Y[, "time"])
    }
    Terms <- attr(m, "terms")
    X <- model.matrix(Terms, m)
    dat <- list(Y = Y, X = X[, -1, drop = FALSE], Xraw = m[, 
        -1, drop = FALSE])
    X <- dat$X
    ##########
    # End of borrowed code    
    ##########                    
                        
                        
                        
    mcmcloop <- mcmcLoop(N=mcmc.control$nits,burnin=mcmc.control$burn,thin=mcmc.control$thin,progressor=mcmcProgressTextBar)                            
                        
    tm <- Y[,1]
    delta <- Y[,2]
    
    mlmod <- flexsurvreg(formula,data=as.data.frame(X),dist="exp")
    estim <- mlmod$opt$par 
    
    betahat <- estim[2:length(estim)]
    omegahat <- estim[1]
    
    SIGMA <- proposalvariance_exponential_nospat(   X=X,
                                                    delta=delta,
                                                    tm=tm,
                                                    betahat=betahat,
                                                    omegahat=omegahat,
                                                    betapriorsd=betapriorsd,
                                                    omegapriorsd=omegapriorsd)                        

    cholsigma <- t(chol(SIGMA))
    
    h <- 1
    
    beta <- betahat
    omega <- omegahat
    d <- length(beta)
    
    oldlogpost <- logposterior_exponential_nospat(tm=tm,delta=delta,X=X,beta=beta,omega=omega,betapriormean=betapriormean,betapriorsd=betapriorsd,omegapriormean=omegapriormean,omegapriorsd=omegapriorsd)
    
    betasamp <- c()
    omegasamp <- c()
    
    tarrec <- rep(NA,mcmc.control$nits)
    
    while(nextStep(mcmcloop)){
    
        newstuff <- c(beta,omega) + h*cholsigma%*%rnorm(d+1)
        
        newlogpost <- logposterior_exponential_nospat(tm=tm,delta=delta,X=X,beta=newstuff[1:d],omega=newstuff[d+1],betapriormean=betapriormean,betapriorsd=betapriorsd,omegapriormean=omegapriormean,omegapriorsd=omegapriorsd)
        
        logfrac <- newlogpost - oldlogpost
        ac <- min(1,exp(logfrac))
        
        if(ac>runif(1)){
            beta <- newstuff[1:d]
            omega <- newstuff[d+1]
            oldlogpost <- newlogpost
        }
        tarrec[iteration(mcmcloop)] <- oldlogpost
        
        h <- exp(log(h) + (1/(iteration(mcmcloop)^0.5))*(ac-0.234))
        
        
        if(is.retain(mcmcloop)){
            betasamp <- rbind(betasamp,as.vector(beta))
            omegasamp <- rbind(omegasamp,as.vector(omega))
        }
    }
    
    #par(mfrow=c(2,1))
    #matplot(betasamp,type="s")
    #matplot(exp(omegasamp),type="s")

    retlist <- list()
    retlist$formula
    retlist$data
    retlist$mcmc.control
    retlist$betapriormean
    retlist$betapriorsd
    retlist$omegapriormean
    retlist$omegapriorsd

    retlist$terms <- Terms
    retlist$mlmod <- mlmod    
    
    retlist$betasamp <- betasamp
    retlist$omegasamp <- omegasamp
    retlist$tarrec <- tarrec
    retlist$lasth <- h
    
    class(retlist) <- c("list","mcmcspatsurv")

    return(retlist)
}

