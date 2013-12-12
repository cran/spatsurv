##' survspat function
##'
##' A function to run a Bayesian analysis on right censored survial data assuming a proportional hazards model with
##' baseline hazard derived from the exponential model.
##'
##' @param formula see ?flexsurvreg 
##' @param data a SpatialPointsDataFrame object
##' @param dist choice of distribution function for baseline hazard options are: "exp"
##' @param covmodel an object of class covmodel, see ?covmodel
##' @param mcmc.control mcmc control parameters, see ?mcmcpars
##' @param priors an object of class Priors, see ?mcmcPriors
##' @return the mcmc output
##' @export


survspat <- function(formula,data,dist,covmodel,mcmc.control,priors){

    start <- Sys.time()

    if(!inherits(data,"SpatialPointsDataFrame")){
        stop("data must be an object of class SpatialPointsDataFrame")
    }
    
    coords <- coordinates(data)
    u <- as.vector(as.matrix(dist(coords)))
    data <- data@data
                        
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
    
    mlmod <- flexsurvreg(formula,data=data,dist=dist)
    estim <- mlmod$opt$par 

    if(dist=="exp"){    
        betahat <- estim[2:length(estim)]
        omegahat <- do.call(paste("transformestimates.",dist,sep=""),args=list(x=estim[1]))
    }
    else if(dist=="weibull"|dist=="gamma"){    
        betahat <- estim[3:length(estim)]
        omegahat <- do.call(paste("transformestimates.",dist,sep=""),args=list(x=estim[1:2]))
    }
    else{
        stop("Unknown dist, must be one of 'exp', 'weibull', or 'gamma'")    
    }
    
    Yhat <- do.call(paste("estimateY.",dist,sep=""),args=list(X=X,betahat=betahat,omegahat=omegahat,tm=tm,delta=delta,u=u,covmodel=covmodel))    
        
    other <- do.call(paste("proposalvariance.",dist,sep=""),args=list(  X=X,
                                                                        delta=delta,
                                                                        tm=tm,
                                                                        betahat=betahat,
                                                                        omegahat=omegahat,
                                                                        Yhat=Yhat,
                                                                        priors=priors,
                                                                        covmodel=covmodel,
                                                                        u=u)) 
    
    gammahat <- other$gammahat
    etahat <- other$etahat                                                                        
    SIGMA <- other$sigma  
                                                                          
    SIGMAINV <- solve(SIGMA)
    cholSIGMA <- Matrix(t(chol(SIGMA)))
    
    h <- 1
    
    beta <- betahat
    omega <- omegahat
    eta <- etahat 
    gamma <- gammahat   
        
    lenbeta <- length(beta)
    lenomega <- length(omega)
    leneta <- length(eta)
    lengamma <- length(gamma)
    
    print(SIGMA[1:(lenbeta+lenomega),1:(lenbeta+lenomega)])
    
    npars <- lenbeta + lenomega + leneta + lengamma
    
    oldlogpost <- do.call(paste("logposterior.",dist,sep=""),args=list(tm=tm,delta=delta,X=X,beta=beta,omega=omega,eta=eta,gamma=gamma,priors=priors,covmodel=covmodel,u=u))
    
    betasamp <- c()
    omegasamp <- c()
    etasamp <- c()
    Ysamp <- c()
    
    tarrec <- oldlogpost$logpost
    
    while(nextStep(mcmcloop)){
    
        stuff <- c(beta,omega,eta,gamma)
        propmean <- stuff + (h/2)*SIGMA%*%oldlogpost$grad
        newstuff <- propmean + h*cholSIGMA%*%rnorm(npars)
        
        newlogpost <- do.call(paste("logposterior.",dist,sep=""),args=list( tm=tm,
                                                                            delta=delta,
                                                                            X=X,
                                                                            beta=newstuff[1:lenbeta],
                                                                            omega=newstuff[(lenbeta+1):(lenbeta+lenomega)],
                                                                            eta=newstuff[(lenbeta+lenomega+1):(lenbeta+lenomega+leneta)],
                                                                            gamma=newstuff[(lenbeta+lenomega+leneta+1):npars],
                                                                            priors=priors,
                                                                            covmodel=covmodel,
                                                                            u=u))
        revmean <- newstuff +  (h/2)*SIGMA%*%newlogpost$grad      
        
        revdiff <- as.matrix(stuff-revmean)
        forwdiff <- as.matrix(newstuff-propmean)
        
        logfrac <- newlogpost$logpost - oldlogpost$logpost -0.5*t(revdiff)%*%SIGMAINV%*%revdiff + 0.5*t(forwdiff)%*%SIGMAINV%*%forwdiff
        ac <- min(1,exp(logfrac))
        
        if(ac>runif(1)){
            beta <- newstuff[1:lenbeta]
            omega <- newstuff[(lenbeta+1):(lenbeta+lenomega)]
            eta <- newstuff[(lenbeta+lenomega+1):(lenbeta+lenomega+leneta)]
            gamma <- newstuff[(lenbeta+lenomega+leneta+1):npars]
            oldlogpost <- newlogpost
        }
        
        h <- exp(log(h) + (1/(iteration(mcmcloop)^0.5))*(ac-0.574))
        if(iteration(mcmcloop)%%100==0){
            cat("\n","h =",h,"\n")
        }
        
        
        if(is.retain(mcmcloop)){
            betasamp <- rbind(betasamp,as.vector(beta))
            omegasamp <- rbind(omegasamp,as.vector(omega))
            etasamp <- rbind(etasamp,as.vector(eta))
            Ysamp <- rbind(Ysamp,as.vector(oldlogpost$Y))
            tarrec <- c(tarrec,oldlogpost$logpost)
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
    retlist$etasamp <- etasamp
    retlist$Ysamp <- Ysamp
    
    retlist$tarrec <- tarrec
    retlist$lasth <- h
    
    retlist$time.taken <- Sys.time() - start
    
    cat("Time taken:",retlist$time.taken,"\n")
    
    class(retlist) <- c("list","mcmcspatsurv")

    return(retlist)
}



##' transformestimates.exp function
##'
##' A function to transform estimates of the parameters of the exponential baseline hazard function, so they are commensurate with R's inbuilt density functions.
##'
##' @param x a vector of paramters 
##' @return the transformed parameters. For the exponential model this is just the identity.
##' @export

transformestimates.exp <- function(x){
    return(x)
}

##' transformestimates.weibull function
##'
##' A function to transform estimates of the parameters of the weibull baseline hazard function, so they are commensurate with R's inbuilt density functions.
##'
##' @param x a vector of paramters
##' @return the transformed parameters. For the weibull mode, this transforms the output from the MALA algorithm so it can be interpreted as 'shape' 'scale', see ?dweibull
##' @export

transformestimates.weibull <- function(x){
    a <- exp(x[1]) # shape
    b <- exp(x[2]) # scale
    alpha <- a
    lambda <- (1/b)^a

    ans <- c(logalpha=log(alpha),loglambda=log(lambda))    
    
    return(ans)
}