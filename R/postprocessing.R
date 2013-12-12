##' print.mcmcspatsurv function
##'
##' A function to print summary tables from an MCMC run 
##'
##' @method print mcmcspatsurv
##' @param x an object inheriting class mcmcspatsurv 
##' @param probs vector of quantiles to return
##' @param digits see help file for format
##' @param scientific see help file for format
##' @param ... additional arguments 
##' @return summary tables
##' @export

print.mcmcspatsurv <- function(x,probs=c(0.5,0.025,0.975),digits = 3, scientific = -3,...){

    quant <- quantile(x,probs)

    cat("\n")
    cat("Fixed Effects:\n")    
    print(quant$betaquant,digits=digits,scientific=scientific)
    cat("\n")
    
    cat("Baseline Hazard Parameters:\n")
    print(quant$omegaquant,digits=digits,scientific=scientific)
    cat("\n")

    cat("Spatial Covariance Parameters:\n")   
    print(quant$etaquant,digits=digits,scientific=scientific)
    cat("\n")
    
    cat("Running Time:\n")
    print(x$time.taken)
} 

##' quantile.mcmcspatsurv function
##'
##' A function to extract quantiles of the parameters from an mcmc run
##'
##' @method quantile mcmcspatsurv
##' @param x an object of class mcmc spatsurv 
##' @param probs vector of probabilities
##' @param ... other arguments to be passed to the function
##' @return ...
##' @export

quantile.mcmcspatsurv <- function(x,probs,...){
    m1 <- t(apply(x$betasamp,2,quantile,probs=probs))
    m2 <- t(apply(x$omegasamp,2,quantile,probs=probs))
    m3 <- t(apply(x$etasamp,2,quantile,probs=probs))
    return(list(betaquant=m1,omegaquant=m2,etaquant=m3))
}



##' summary.mcmcspatsurv function
##'
##' A function to print summary tables from an MCMC run 
##'
##' @method summary mcmcspatsurv
##' @param object an object inheriting class mcmcspatsurv 
##' @param probs vector of quantiles to return
##' @param ... additional arguments 
##' @return summary tables
##' @export

summary.mcmcspatsurv <- function(object,probs=c(0.5,0.025,0.975),...){
    quant <- quantile(object,probs)
    return(rbind(quant$betaquant,quant$omegaquant,quant$etaquant))
} 


##' frailtylag1 function
##'
##' A function to produce and return the lag 1 autocorrelation for each of the spatially correlated frailty chains
##'
##' @param object an object of class mcmcspatsurv 
##' @param ... other arguments to be passed to the plot function 
##' @return the lag 1 autocorrelation for each of the spatially correlated frailty chains 
##' @export

frailtylag1 <- function(object,...){
    lag1acf <- apply(object$Ysamp,2,function(x){acf(x,plot=FALSE)$acf[2]})
    plot(lag1acf,xlab="Frailty Index",ylab="Lag 1 Autocorrelation",ylim=c(-1,1),...)
    return(lag1acf)
}




##' spatialpars function
##'
##' A function to return the mcmc chains for the spatial covariance function parameters
##'
##' @param x an object of class mcmcspatsurv
##' @return the mcmc chains
##' @export

spatialpars <- function(x){
    return(x$etasamp)
}



##' hazardpars function
##'
##' A function to return the mcmc chains for the hazard function parameters
##'
##' @param x an object of class mcmcspatsurv
##' @return the mcmc chains
##' @export

hazardpars <- function(x){
    return(x$omegasamp)
}



##' fixedpars function
##'
##' A function to return the mcmc chains for the covariate effects
##'
##' @param x an object of class mcmcspatsurv
##' @return the mcmc chains
##' @export

fixedpars <- function(x){
    return(x$etasamp)
}



##' randompars function
##'
##' A function to return the mcmc chains for the spatially correlated frailties
##'
##' @param x an object of class mcmcspatsurv
##' @return the mcmc chains
##' @export

randompars <- function(x){
    return(x$etasamp)
}



##' baselinehazard function
##'
##' A function to 
##'
##' @param x X 
##' @param t X 
##' @param n X
##' @param probs X
##' @param plot X 
##' @return ...
##' @export

baselinehazard <- function(x,t=NULL,n=100,probs=c(0.025,0.5,0.975),plot=TRUE){
    # extract and transform onto appropriate scale    
    transfun <- get(paste("transformestimates.",x$dist,sep=""))
    if(ncol(x$omegasamp)==1){    
        omegasamp <- matrix(apply(x$omegasamp,1,transfun))
    }
    else{
        omegasamp <- t(apply(x$omegasamp,1,transfun))
    }    
    
    if(is.null(t)){
        t <- seq(0,max(x$mlmod$data$Y[,"time"]),length.out=n)
    }
    
    fun <- function(pars){
        f <- get(paste("basehaz.",x$dist,sep=""))(pars)
        return(f(t))
    } 
    
    samp <- t(apply(omegasamp,1,fun))   

    toreturn <- t(apply(samp,2,quantile,probs=probs))
    
    rownames(toreturn) <- t 
    
    if(plot){
        if(length(probs)==3){
            matplot(toreturn,type="l",col=c("purple","black","blue"),lty=c("dashed","solid","dashed"),xlab="time",ylab="Baseline Hazard")
            legend("topright",lty=c("dashed","solid","dashed"),col=rev(c("purple","black","blue")),legend=rev(probs))
        }
        else{
            matplot(toreturn,type="l",xlab="time",ylab="Baseline Hazard")
        }
    }    
    
    return(toreturn)
}


##' basehaz.exp function
##'
##' A function to 
##'
##' @param pars X 
##' @return ...
##' @export

basehaz.exp <- function(pars){
    f <- function(t){
        return(pars)
    }
    return(f)
}



##' basehaz.weibull function
##'
##' A function to 
##'
##' @param pars X 
##' @return ...
##' @export

basehaz.weibull <- function(pars){
    alpha <- pars[1]
    lambda <- pars[2]
    f <- function(t){
        return(lambda*alpha*t^(alpha-1))
    }
    return(f)
}



##' hazard_exp function
##'
##' A function to compute the hazard function for an individual where the baseline hazard comes from an exponential survival model
##'
##' @param inputs inputs for the function including the model matrix, frailties, fixed effects and the parameters of the baseline hazard derived from this model
##' @return the hazard function for the individual
##' @export

hazard_exp <- function(inputs){
    X <- inputs$X
    Y <- inputs$Y
    beta <- inputs$beta
    
    Xbeta <- sum(X*beta)
    expXbeta <- exp(Xbeta)
    expXbeta_plus_Y <- expXbeta*exp(Y)
    
    rate <- inputs$omega
    
    f <- function(t){
        return(expXbeta_plus_Y*rate)
    }
    return(f)      
}




##' survival_exp function
##'
##' A function to compute the survival function for an individual where the baseline hazard comes from an exponential survival model
##'
##' @param inputs inputs for the function including the model matrix, frailties, fixed effects and the parameters of the baseline hazard derived from this model
##' @return the survival function for the individual
##' @export
survival_exp <- function(inputs){
    X <- inputs$X
    Y <- inputs$Y
    beta <- inputs$beta
    
    Xbeta <- sum(X*beta)
    expXbeta <- exp(Xbeta)
    expXbeta_plus_Y <- expXbeta*exp(Y)
    
    rate <- inputs$omega
    
    f <- function(t){
        return(exp(-expXbeta_plus_Y*rate*t))
    }
    return(f)      
}



##' density_exp function
##'
##' A function to compute the density function for an individual where the baseline hazard comes from an exponential survival model
##'
##' @param inputs inputs for the function including the model matrix, frailties, fixed effects and the parameters of the baseline hazard derived from this model
##' @return the density function for the individual
##' @export

density_exp <- function(inputs){
    X <- inputs$X
    Y <- inputs$Y
    beta <- inputs$beta
    
    Xbeta <- sum(X*beta)
    expXbeta <- exp(Xbeta)
    expXbeta_plus_Y <- expXbeta*exp(Y)
    
    rate <- inputs$omega
    
    f <- function(t){
        return(expXbeta_plus_Y*rate*exp(-expXbeta_plus_Y*rate*t))
    }
    return(f)      
}



##' densityquantile_exp function
##'
##' A function to compute quantiles of the density function for an individual where the baseline hazard comes from an exponential survival model
##'
##' @param inputs inputs for the function including the model matrix, frailties, fixed effects and the parameters of the baseline hazard derived from this model
##' @return quantiles of the density function for the individual
##' @export

densityquantile_exp <- function(inputs){
    X <- inputs$X
    Y <- inputs$Y
    beta <- inputs$beta
    
    Xbeta <- sum(X*beta)
    expXbeta <- exp(Xbeta)
    expXbeta_plus_Y <- expXbeta*exp(Y)
    
    rate <- inputs$omega
    
    f <- function(probs){
        a <- expXbeta_plus_Y*rate
        return((-1/a)*log(1-probs))
    }
    return(f)    
}




##' hazard_weibull function
##'
##' A function to compute the hazard function for an individual where the baseline hazard comes from a Weibull survival model
##'
##' @param inputs inputs for the function including the model matrix, frailties, fixed effects and the parameters of the baseline hazard derived from this model
##' @return the hazard function for the individual
##' @export

hazard_weibull <- function(inputs){
    X <- inputs$X
    Y <- inputs$Y
    beta <- inputs$beta
    
    Xbeta <- sum(X*beta)
    expXbeta <- exp(Xbeta)
    expXbeta_plus_Y <- expXbeta*exp(Y)
    
    alpha <- inputs$omega[1]
    lambda <- inputs$omega[2] 
    
    f <- function(t){
        return(expXbeta_plus_Y*lambda*alpha*t^(alpha-1))
    }
    return(f)    
}



##' survival_weibull function
##'
##' A function to compute the survival function for an individual where the baseline hazard comes from a Weibull survival model
##'
##' @param inputs inputs for the function including the model matrix, frailties, fixed effects and the parameters of the baseline hazard derived from this model
##' @return the survival function for the individual
##' @export

survival_weibull <- function(inputs){
    X <- inputs$X
    Y <- inputs$Y
    beta <- inputs$beta
    
    Xbeta <- sum(X*beta)
    expXbeta <- exp(Xbeta)
    expXbeta_plus_Y <- expXbeta*exp(Y)
    
    alpha <- inputs$omega[1]
    lambda <- inputs$omega[2]   
    
    f <- function(t){
        return(exp(-expXbeta_plus_Y*lambda*t^alpha))
    }
    return(f)    
}



##' density_weibull function
##'
##' A function to compute the density function for an individual where the baseline hazard comes from a Weibull survival model
##'
##' @param inputs inputs for the function including the model matrix, frailties, fixed effects and the parameters of the baseline hazard derived from this model
##' @return the density function for the individual
##' @export

density_weibull <- function(inputs){
    X <- inputs$X
    Y <- inputs$Y
    beta <- inputs$beta
    
    Xbeta <- sum(X*beta)
    expXbeta <- exp(Xbeta)
    expXbeta_plus_Y <- expXbeta*exp(Y)
    
    alpha <- inputs$omega[1]
    lambda <- inputs$omega[2]
    
    f <- function(t){
        return(expXbeta_plus_Y*lambda*alpha*t^(alpha-1)*exp(-expXbeta_plus_Y*lambda*t^alpha))
    }
    return(f)    
}



##' densityquantile_weibull function
##'
##' A function to compute quantiles of the density function for an individual where the baseline hazard comes from a Weibull survival model
##'
##' @param inputs inputs for the function including the model matrix, frailties, fixed effects and the parameters of the baseline hazard derived from this model
##' @return quantiles of the density function for the individual
##' @export

densityquantile_weibull <- function(inputs){
    X <- inputs$X
    Y <- inputs$Y
    beta <- inputs$beta
    
    Xbeta <- sum(X*beta)
    expXbeta <- exp(Xbeta)
    expXbeta_plus_Y <- expXbeta*exp(Y)
    
    alpha <- inputs$omega[1]
    lambda <- inputs$omega[2]
    
    f <- function(probs){
        a <- expXbeta_plus_Y*lambda
        return(((-1/a)*log(1-probs))^(1/alpha))
    }
    return(f)    
}


##' predict.mcmcspatsurv function
##'
##' A function to 
##'
##' @method predict mcmcspatsurv
##' @param object an object of class mcmcspatsurv
##' @param type can be "densityquantile","density", "hazard" or "survival". Default is "densityquantile".
##' @param newdata X
##' @param t X
##' @param n X
##' @param probs X
##' @param ask X
##' @param ... other arguments 
##' @return ...
##' @export

predict.mcmcspatsurv <- function(object,type="densityquantile",newdata=NULL,t=NULL,n=110,probs=c(0.025,0.5,0.975),ask=TRUE,...){

    if(is.null(newdata)){
        newdata <- object$mlmod$data$X
    }
    nobs <- nrow(newdata)
    
    if(is.null(t)){
        t <- seq(0,max(object$mlmod$data$Y[,"time"]),length.out=n)
    }

    predictmat <- NULL    
    if(type=="densityquantile"){
        t <- probs
        n <- length(t)
        pb <- txtProgressBar(min = 1, max = nobs)
        predictmat <- matrix(NA,nobs,length(probs))
    }
    
    nits <- nrow(object$Ysamp)

        
    if(type!="densityquantile"){
        par(ask=ask)
    }
        
    
    for(i in 1:nobs){
        dat <- matrix(NA,nits,n)
        inputs <- list()
        inputs$X <- newdata[i,]
        for (j in 1:nits){
            inputs$Y <- object$Ysamp[j,i]
            inputs$beta <- object$betasamp[j,]
            inputs$omega <- object$omegasamp[j,]
            fun <- get(paste(type,"_",object$dist,sep=""))(inputs)
            dat[j,] <- fun(t)      
        }

        if(type!="densityquantile"){            
            toplot <- t(apply(dat,2,quantile,probs=probs))        
            matplot(t,toplot,type="l",col=c("purple","black","blue"),lty=c("dashed","solid","dashed"),xlab="time",ylab=type,main=paste("Individual",i))
            legend("topright",lty=c("dashed","solid","dashed"),col=rev(c("purple","black","blue")),legend=rev(probs))
        }
        else{        
            setTxtProgressBar(pb,i)
        }
        
        if(!is.null(predictmat)){
            predictmat[i,] <- colMeans(dat)
        }       
    }
    
    if(type=="densityquantile"){
        close(pb)
    }   
   
    return(predictmat)
    
}


## plot.mcmcspatsurv function
##
## A function to produce diagnostic plots for objects of class mcmcspatsurv
##
## @method plot mcmcspatsurv
## @param x an object of class mcmcspatsurv
## @param n number of time points to consider
## @param pr optional predictions, if they've already been computed, must be type="densityquantile", see ?predict.mcmcspatsurv. If NULL, the predicted median survival time is used.
## @param alpha significance level, default is 0.05. 
## @param ... other arguments  
## @return produces some diagnostic plots (currently only one diagnostic plot...)
## @export
#
#plot.mcmcspatsurv <- function(x,n=1000,pr=NULL,alpha=0.05,...){
#    survdat <- getsurvdata(x)
#    times <- survdat[,1]
#    cens <- survdat[,2]
#    maxt <- max(times)
#    tm <- seq(min(times),maxt,length.out=n)
#    
#    if(is.null(pr)){
#        cat("To save time you can use the pr argument to feed in precomputed predictions using the predict function, see ?predict.mcmcspatsurv.\n")
#        pr <- predict(x,type="densityquantile",probs=0.5)
#    }
#    
#    p <- rep(0,n)
#    np <- rep(0,n)
#    o <- rep(0,n)
#    no <- rep(0,n)
#    for(i in 1:n){
#        temptimes <- times
#        temptimes[times<=tm[i] & cens==0] <- NA
#        p[i] <- sum(pr>tm[i])
#        np[i] <- length(times)
#        o[i] <- sum(temptimes>tm[i] & cens==1,na.rm=TRUE)
#        no[i] <- sum(!is.na(temptimes))
#    }
#    
#    rr <- (p/np)/(o/no)
#    logrr <- log(rr)
#    
#    selogrr <- sqrt(1/p-1/np+1/o-1/no)
#    z <- qnorm(1-alpha/2)
#    lower <- exp(logrr - z*selogrr)
#    upper <- exp(logrr + z*selogrr)
#    
#    lower[is.na(lower)|is.infinite(lower)] <- NA
#    upper[is.na(upper)|is.infinite(upper)] <- NA
#    
#    tmax <- tm[min(which(is.na(lower))[1],which(is.na(upper))[1])]
#    
#    plot(NULL,xlim=c(0,tmax),ylim=range(c(lower,upper),na.rm=TRUE),xlab="Time",ylab="propn. predicted to survive / propn. observed to survive")
#    lines(tm,lower,col="red",lty="dashed")
#    lines(tm,upper,col="red",lty="dashed")
#    lines(tm,rr)
#    abline(h=1,col="blue")
#    return(list(pr=pr,tm=tm,rr=rr,selogrr=selogrr,lower=lower,upper=upper))
#}



priorposterior <- function(x,breaks=30,ylab="Density",main="",ask=TRUE,...){
    nbeta <- ncol(x$betasamp)
    nomega <- ncol(x$omegasamp)
    neta <- ncol(x$etasamp)
    
    par(ask=ask)
    
    for(i in 1:nbeta){
        h <- hist(x$betasamp[,i],ylab=ylab,xlab=colnames(x$betasamp)[i],breaks=breaks,freq=FALSE,main=main,...)
        xrg <- range(h$breaks)
        r <- seq(xrg[1],xrg[2],length.out=1000)
        if(length(x$priors$betaprior$mean)==1){
            lines(r,dnorm(r,mean=x$priors$betaprior$mean,sd=x$priors$betaprior$sd),col="red",lwd=2)
        }
        else{
            lines(r,dnorm(r,mean=x$priors$betaprior$mean[i],sd=x$priors$betaprior$sd[i]),col="red",lwd=2)
        }
    }
    
    
    transfun <- get(paste("transformestimates.",x$dist,sep=""))
    samp <- t(apply(x$omegasamp,1,transfun))
    samp <- t(apply(samp,1,x$omegatrans))
    for(i in 1:nomega){        
        h <- hist(samp[,i],ylab=ylab,xlab=paste("Transformed",colnames(x$omegasamp)[i]),breaks=breaks,freq=FALSE,main=main,...)
        xrg <- range(h$breaks)
        r <- seq(xrg[1],xrg[2],length.out=1000)
        if(length(x$priors$betaprior$mean)==1){
            lines(r,dnorm(r,mean=x$priors$omegaprior$mean,sd=x$priors$omegaprior$sd),col="red",lwd=2)
        }
        else{
            lines(r,dnorm(r,mean=x$priors$omegaprior$mean[i],sd=x$priors$omegaprior$sd[i]),col="red",lwd=2)
        }
    }
    
    for(i in 1:neta){
        samp <- x$cov.model$trans[[i]](x$etasamp[,i])
        h <- hist(samp,ylab=ylab,xlab=paste("Transformed",colnames(x$etasamp)[i]),breaks=breaks,freq=FALSE,main=main,...)
        xrg <- range(h$breaks)
        r <- seq(xrg[1],xrg[2],length.out=1000)
        if(length(x$priors$etaprior$mean)==1){
            lines(r,dnorm(r,mean=x$priors$etaprior$mean,sd=x$priors$etaprior$sd),col="red",lwd=2)
        }
        else{
            lines(r,dnorm(r,mean=x$priors$etaprior$mean[i],sd=x$priors$etaprior$sd[i]),col="red",lwd=2)
        }
    }
}



##' posteriorcov function
##'
##' A function to produce a plot of the posterior covariance function.
##'
##' @param x an object of class mcmcspatsurv 
##' @param probs vector of probabilities to be fed to quantile function  
##' @param rmax  maximum distance in space to compute this distance up to
##' @param n the number of points at which to evaluate the posterior covariance.
##' @param ... other arguments to be passed to matplot function 
##' @return produces a plot of the posterior spatial covariance function.
##' @export

posteriorcov <- function(x,probs=c(0.025,0.5,0.975),rmax=NULL,n=100,...){
    nr <- nrow(x$etasamp)
    nc <- ncol(x$etasamp)
    pars <- matrix(NA,nr,nc)
    for(i in 1:nc){
        pars[,i] <- x$cov.model$trans[[i]](x$etasamp[,i]) # transform (e.g. to log-scal)
    }
    rmaxx <- 0.25*sum(apply(bbox(x$data),1,diff))/2 # approx 1/4 of mean length of observation window
    if(!is.null(rmax)){
        rmaxx <- rmax
    }
    
    r <- seq(0,rmaxx,length.out=n)
    covs <- t(apply(pars,1,function(x){x$cov.model$eval(r,pars=x)}))
    qts <- t(apply(covs,2,quantile,probs=probs))
    
    matplot(r,qts,type="l",xlab="Distance",ylab="Covariance",...)
}