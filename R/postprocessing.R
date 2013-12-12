##' print.mcmcspatsurv function
##'
##' A function to print summary tables from an MCMC run 
##'
##' @method print mcmcspatsurv
##' @param x an object inheriting class mcmcspatsurv 
##' @param probs vector of quantiles to return
##' @param digits see help file ?format
##' @param scientific see help file ?format
##' @param ... additional arguments, not used here 
##' @return prints summary tables to the console
##' @seealso \link{quantile.mcmcspatsurv}, \link{summary.mcmcspatsurv}, \link{vcov.mcmcspatsurv}, 
##' \link{frailtylag1}, \link{spatialpars}, \link{hazardpars}, \link{fixedpars}, \link{randompars},
##' \link{baselinehazard}, \link{predict.mcmcspatsurv}, \link{priorposterior}, \link{posteriorcov}, \link{MCE},
##' \link{hazardexceedance}
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
    
    cat("Deviance Information Criterion: ",x$DIC,"\n")
    cat("\n")
    
    cat("MCMC Details:\n")
    m <- matrix(c(x$mcmc.control$nits,x$mcmc.control$burn,x$mcmc.control$thin),3,1)
    rownames(m) <- c("Number of Iterations","Burnin length","Thining")
    colnames(m) <- ""
    print(m,digits=digits,scientific=scientific)
    cat("\n")
    
    cat("Running Time:\n")
    print(x$time.taken)
} 

##' quantile.mcmcspatsurv function
##'
##' A function to extract quantiles of the parameters from an mcmc run
##'
##' @method quantile mcmcspatsurv
##' @param x an object inheriting class mcmcspatsurv 
##' @param probs vector of probabilities
##' @param ... other arguments to be passed to the function, not used here
##' @return quantiles of model parameters
##' @seealso \link{print.mcmcspatsurv}, \link{summary.mcmcspatsurv}, \link{vcov.mcmcspatsurv}, 
##' \link{frailtylag1}, \link{spatialpars}, \link{hazardpars}, \link{fixedpars}, \link{randompars},
##' \link{baselinehazard}, \link{predict.mcmcspatsurv}, \link{priorposterior}, \link{posteriorcov}, \link{MCE},
##' \link{hazardexceedance}
##' @export

quantile.mcmcspatsurv <- function(x,probs=c(0.025,0.5,0.975),...){
    m1 <- t(apply(x$betasamp,2,quantile,probs=probs))
    m2 <- t(apply(x$omegasamp,2,quantile,probs=probs))
    m3 <- t(apply(x$etasamp,2,quantile,probs=probs))
    return(list(betaquant=m1,omegaquant=m2,etaquant=m3))
}



##' summary.mcmcspatsurv function
##'
##' A function to return summary tables from an MCMC run 
##'
##' @method summary mcmcspatsurv
##' @param object an object inheriting class mcmcspatsurv 
##' @param probs vector of quantiles to return
##' @param ... additional arguments 
##' @return summary tables to the console
##' @seealso \link{print.mcmcspatsurv}, \link{quantile.mcmcspatsurv}, \link{vcov.mcmcspatsurv}, 
##' \link{frailtylag1}, \link{spatialpars}, \link{hazardpars}, \link{fixedpars}, \link{randompars},
##' \link{baselinehazard}, \link{predict.mcmcspatsurv}, \link{priorposterior}, \link{posteriorcov}, \link{MCE},
##' \link{hazardexceedance}
##' @export

summary.mcmcspatsurv <- function(object,probs=c(0.5,0.025,0.975),...){
    quant <- quantile(object,probs)
    return(rbind(quant$betaquant,quant$omegaquant,quant$etaquant))
} 



##' vcov.mcmcspatsurv function
##'
##' A function to return the variance covariance matrix of the parameters beta, omega and eta
##'
##' @method vcov mcmcspatsurv
##' @param object an object inheriting class mcmcspatsurv 
##' @param ... other arguments, not used here
##' @return the variance covariance matrix of the parameters beta, omega and eta
##' @seealso \link{print.mcmcspatsurv}, \link{quantile.mcmcspatsurv}, \link{summary.mcmcspatsurv}, 
##' \link{frailtylag1}, \link{spatialpars}, \link{hazardpars}, \link{fixedpars}, \link{randompars},
##' \link{baselinehazard}, \link{predict.mcmcspatsurv}, \link{priorposterior}, \link{posteriorcov}, \link{MCE},
##' \link{hazardexceedance}
##' @export

vcov.mcmcspatsurv <- function(object,...){
    return(cov(cbind(object$betasamp,object$omegasamp,object$etasamp)))
} 



##' frailtylag1 function
##'
##' A function to produce a plot of, and return, the lag 1 autocorrelation for each of the spatially correlated frailty chains
##'
##' @param object an object inheriting class mcmcspatsurv 
##' @param plot logical whether to plot the result, default is TRUE
##' @param ... other arguments to be passed to the plot function 
##' @return the lag 1 autocorrelation for each of the spatially correlated frailty chains
##' @seealso \link{print.mcmcspatsurv}, \link{quantile.mcmcspatsurv}, \link{summary.mcmcspatsurv}, \link{vcov.mcmcspatsurv}, 
##' \link{spatialpars}, \link{hazardpars}, \link{fixedpars}, \link{randompars},
##' \link{baselinehazard}, \link{predict.mcmcspatsurv}, \link{priorposterior}, \link{posteriorcov}, \link{MCE},
##' \link{hazardexceedance} 
##' @export

frailtylag1 <- function(object,plot=TRUE,...){
    lag1acf <- apply(object$Ysamp,2,function(x){acf(x,plot=FALSE)$acf[2]})
    if(plot){
        plot(lag1acf,xlab="Frailty Index",ylab="Lag 1 Autocorrelation",ylim=c(-1,1),...)
    }
    return(lag1acf)
}




##' spatialpars function
##'
##' A function to return the mcmc chains for the spatial covariance function parameters
##'
##' @param x an object of class mcmcspatsurv
##' @return the eta mcmc chains
##' @seealso \link{print.mcmcspatsurv}, \link{quantile.mcmcspatsurv}, \link{summary.mcmcspatsurv}, \link{vcov.mcmcspatsurv}, 
##' \link{frailtylag1}, \link{hazardpars}, \link{fixedpars}, \link{randompars},
##' \link{baselinehazard}, \link{predict.mcmcspatsurv}, \link{priorposterior}, \link{posteriorcov}, \link{MCE},
##' \link{hazardexceedance}
##' @export

spatialpars <- function(x){
    return(x$etasamp)
}



##' hazardpars function
##'
##' A function to return the mcmc chains for the hazard function parameters
##'
##' @param x an object of class mcmcspatsurv
##' @return the omega mcmc chains
##' @seealso \link{print.mcmcspatsurv}, \link{quantile.mcmcspatsurv}, \link{summary.mcmcspatsurv}, \link{vcov.mcmcspatsurv}, 
##' \link{frailtylag1}, \link{spatialpars}, \link{fixedpars}, \link{randompars},
##' \link{baselinehazard}, \link{predict.mcmcspatsurv}, \link{priorposterior}, \link{posteriorcov}, \link{MCE},
##' \link{hazardexceedance}
##' @export

hazardpars <- function(x){
    return(x$omegasamp)
}



##' fixedpars function
##'
##' A function to return the mcmc chains for the covariate effects
##'
##' @param x an object of class mcmcspatsurv
##' @return the beta mcmc chains
##' @seealso \link{print.mcmcspatsurv}, \link{quantile.mcmcspatsurv}, \link{summary.mcmcspatsurv}, \link{vcov.mcmcspatsurv}, 
##' \link{frailtylag1}, \link{spatialpars}, \link{hazardpars}, \link{randompars},
##' \link{baselinehazard}, \link{predict.mcmcspatsurv}, \link{priorposterior}, \link{posteriorcov}, \link{MCE},
##' \link{hazardexceedance}
##' @export

fixedpars <- function(x){
    return(x$betasamp)
}



##' randompars function
##'
##' A function to return the mcmc chains for the spatially correlated frailties
##'
##' @param x an object of class mcmcspatsurv
##' @return the Y mcmc chains
##' @seealso \link{print.mcmcspatsurv}, \link{quantile.mcmcspatsurv}, \link{summary.mcmcspatsurv}, \link{vcov.mcmcspatsurv}, 
##' \link{frailtylag1}, \link{spatialpars}, \link{hazardpars}, \link{fixedpars},
##' \link{baselinehazard}, \link{predict.mcmcspatsurv}, \link{priorposterior}, \link{posteriorcov}, \link{MCE},
##' \link{hazardexceedance}
##' @export

randompars <- function(x){
    return(x$Ysamp)
}



##' baselinehazard function
##'
##' A function to compute quantiles of the posterior baseline hazard or cumulative baseline hazard.
##'
##' @param x an object inheriting class mcmcspatsurv 
##' @param t optional vector of times at which to compute the quantiles, Defult is NULL, in which case a uniformly spaced vector of length n from 0 to the maximum time is used 
##' @param n the number of points at which to compute the quantiles if t is NULL
##' @param probs vector of probabilities 
##' @param cumulative logical, whether to return the baseline hazard (default i.e. FALSE) or cumulative baseline hazard 
##' @param plot whether to plot the result
##' @return the vector of times and quantiles of the baseline or cumulative baseline hazard at those times
##' @seealso \link{print.mcmcspatsurv}, \link{quantile.mcmcspatsurv}, \link{summary.mcmcspatsurv}, \link{vcov.mcmcspatsurv}, 
##' \link{frailtylag1}, \link{spatialpars}, \link{hazardpars}, \link{fixedpars}, \link{randompars},
##' \link{predict.mcmcspatsurv}, \link{priorposterior}, \link{posteriorcov}, \link{MCE},
##' \link{hazardexceedance}
##' @export

baselinehazard <- function(x,t=NULL,n=100,probs=c(0.025,0.5,0.975),cumulative=FALSE,plot=TRUE){

    omegasamp <- x$omegasamp   
    
    if(is.null(t)){
        if(x$censoringtype=="left" | x$censoringtype=="right"){
            t <- seq(0,max(x$survivaldata[,"time"],na.rm=TRUE),length.out=n)
        }
        else{
            t <- seq(0,max(c(x$survivaldata[,"time1"],x$survivaldata[,"time1"]),na.rm=TRUE),length.out=n)
        }
    }
    
    fun <- function(pars){
        f <- basehazard(x$dist)(pars)
        return(f(t))
    }
    
    YLAB <- "Baseline Hazard"
    if(cumulative){ 
        fun <- function(pars){
            f <- cumbasehazard(x$dist)(pars)
            return(f(t))
        }
        YLAB <- "Cumulative Baseline Hazard"
    }    
    samp <- t(apply(omegasamp,1,fun))   

    toreturn <- t(apply(samp,2,quantile,probs=probs))
    
    rownames(toreturn) <- t 
    
    if(plot){
        if(length(probs)==3){
            matplot(t,toreturn,type="l",col=c("purple","black","blue"),lty=c("dashed","solid","dashed"),xlab="time",ylab=YLAB)
            legend("topright",lty=c("dashed","solid","dashed"),col=rev(c("purple","black","blue")),legend=rev(probs))
        }
        else{
            matplot(t,toreturn,type="l",xlab="time",ylab="Baseline Hazard")
        }
    }    
    
    return(list(t=t,qts=toreturn))
}






##' hazard_PP function
##'
##' A function to compute an individual's hazard function.
##'
##' @param inputs inputs for the function including the model matrix, frailties, fixed effects and the parameters of the baseline hazard derived from this model
##' @return the hazard function for the individual
##' @export

hazard_PP <- function(inputs){
    X <- inputs$X
    Y <- inputs$Y
    beta <- inputs$beta
    
    Xbeta <- sum(X*beta)
    expXbeta <- exp(Xbeta)
    expXbeta_plus_Y <- expXbeta*exp(Y)
    
    f <- function(t){
        h <- basehazard(inputs$dist)(inputs$omega)
        return(expXbeta_plus_Y*h(t))
    }
    return(f)      
}




##' survival_PP function
##'
##' A function to compute an individual's survival function 
##'
##' @param inputs inputs for the function including the model matrix, frailties, fixed effects and the parameters of the baseline hazard derived from this model
##' @return the survival function for the individual
##' @export
survival_PP <- function(inputs){
    X <- inputs$X
    Y <- inputs$Y
    beta <- inputs$beta
    
    Xbeta <- sum(X*beta)
    expXbeta <- exp(Xbeta)
    expXbeta_plus_Y <- expXbeta*exp(Y)
    
    f <- function(t){
        H <- cumbasehazard(inputs$dist)(inputs$omega)
        return(exp(-expXbeta_plus_Y*H(t)))
    }
    return(f)      
}



##' density_PP function
##'
##' A function to compute an individual's density function 
##'
##' @param inputs inputs for the function including the model matrix, frailties, fixed effects and the parameters of the baseline hazard derived from this model
##' @return the density function for the individual
##' @export

density_PP <- function(inputs){
    X <- inputs$X
    Y <- inputs$Y
    beta <- inputs$beta
    
    Xbeta <- sum(X*beta)
    expXbeta <- exp(Xbeta)
    expXbeta_plus_Y <- expXbeta*exp(Y)
    
    f <- function(t){
        h <- basehazard(inputs$dist)(inputs$omega)    
        H <- cumbasehazard(inputs$dist)(inputs$omega)
        return(expXbeta_plus_Y*h(t)*exp(-expXbeta_plus_Y*H(t)))
    }
    return(f)      
}



##' Et_PP function
##'
##' A function to compute an individual's approximate expected survival time using numerical integration. Note this appears to be unstable; the
##' function is based on R's integrate function. Not intended for general use (yet!).  
##'
##' @param inputs inputs for the function including the model matrix, frailties, fixed effects and the parameters of the baseline hazard derived from this model
##' @return the expected survival time for the individual, obtained by numerical integration of the density function.
##' @export

Et_PP <- function(inputs){
    f <- density_PP(inputs)
    expect <- function(t){
        return(t*f(t))
    }
    int <- try(integrate(expect,0,Inf)$value,silent=TRUE)
    if(inherits(int,"try-error")){
        return(NA)
    }
    else{
        return(int)
    } 
}



##' densityquantile_PP function
##'
##' A function to compute quantiles of the density function 
##'
##' @param inputs inputs for the function including the model matrix, frailties, fixed effects and the parameters of the baseline hazard derived from this model
##' @return quantiles of the density function for the individual
##' @export

densityquantile_PP <- function(inputs){
    X <- inputs$X
    Y <- inputs$Y
    beta <- inputs$beta
    
    Xbeta <- sum(X*beta)
    expXbeta <- exp(Xbeta)
    expXbeta_plus_Y <- expXbeta*exp(Y)
    
    dq <- densityquantile(inputs$dist)(inputs$omega,other=list(expXbetaplusY=expXbeta_plus_Y)) 
    
    return(dq)    
}





##' predict.mcmcspatsurv function
##'
##' A function to produce predictions from MCMC output. These could include quantiles of the individual density, survival or 
##' hazard functions or quantiles of the density function (if available analytically).
##'
##' @method predict mcmcspatsurv
##' @param object an object of class mcmcspatsurv
##' @param type can be "density", "hazard", "survival" or "densityquantile". Default is "density". Note that "densityquantile" is not always analytically tractable for some choices of baseline hazard function.
##' @param t optional vector of times at which to compute the quantiles, Defult is NULL, in which case a uniformly spaced vector of length n from 0 to the maximum time is used
##' @param n the number of points at which to compute the quantiles if t is NULL
##' @param indx the index number of a particular individual or vector of indices of individuals for which the quantiles should be produced 
##' @param probs vector of probabilities
##' @param plot whether to plot the result
##' @param pause logical whether to pause between plots, the default is TRUE
##' @param ... other arguments, not used here 
##' @return the required predictions
##' @seealso \link{print.mcmcspatsurv}, \link{quantile.mcmcspatsurv}, \link{summary.mcmcspatsurv}, \link{vcov.mcmcspatsurv}, 
##' \link{frailtylag1}, \link{spatialpars}, \link{hazardpars}, \link{fixedpars}, \link{randompars},
##' \link{baselinehazard}, \link{priorposterior}, \link{posteriorcov}, \link{MCE},
##' \link{hazardexceedance}
##' @export

predict.mcmcspatsurv <- function(object,type="density",t=NULL,n=110,indx=NULL,probs=c(0.025,0.5,0.975),plot=TRUE,pause=TRUE,...){

    newdata <- object$X

    if(is.null(indx)){
        indx <- 1:nrow(newdata)
    }    
    
    nobs <- length(indx) # nrow(newdata)
    
    if(is.null(t)){
        if(object$censoringtype=="left" | object$censoringtype=="right"){
            t <- seq(0,max(object$survivaldata[,"time"],na.rm=TRUE),length.out=n)
        }
        else{
            t <- seq(0,max(c(object$survivaldata[,"time1"],object$survivaldata[,"time1"]),na.rm=TRUE),length.out=n)
        }
    }

    predictmat <- NULL    
    if(type=="densityquantile"){
        t <- probs
        n <- length(t)
        pb <- txtProgressBar(min = 0, max = nobs)
        predictmat <- matrix(NA,nobs,length(probs))
    }
    
    if(type=="Et"){
        n <- nobs
        pb <- txtProgressBar(min = 0, max = nobs)
    }
    
    
    
    nits <- nrow(object$Ysamp)
    
    omegasamp <- object$omegasamp    
    dat <- matrix(NA,nits,n)
    for(i in indx){
        if(type!="Et"){
            dat <- matrix(NA,nits,n)
        }
        inputs <- list()
        inputs$X <- newdata[i,]
        inputs$dist <- object$dist
        for (j in 1:nits){
            inputs$Y <- object$Ysamp[j,i]
            inputs$beta <- object$betasamp[j,]
            inputs$omega <- omegasamp[j,]
            fun <- get(paste(type,"_PP",sep=""))(inputs)
            if(type!="Et"){
                dat[j,] <- fun(t)
            }
            else{
                dat[j,i] <- fun
            }      
        }

        if(type!="densityquantile" & type!="Et"){          
            toplot <- t(apply(dat,2,quantile,probs=probs))        
            matplot(t,toplot,type="l",col=c("purple","black","blue"),lty=c("dashed","solid","dashed"),xlab="time",ylab=type,main=paste("Individual",i))
            legend("topright",lty=c("dashed","solid","dashed"),col=rev(c("purple","black","blue")),legend=rev(probs))
            if(pause){
                cat("[press [enter] to continue]")
                ob <- scan(n=1,quiet=TRUE)
            }
        }
        else{        
            setTxtProgressBar(pb,i)
        }
        
        if(!is.null(predictmat)){
            if(type!="Et"){
                predictmat[i,] <- colMeans(dat)
            }
        }     
    }
    
    if(type=="densityquantile" & type=="Et"){
        close(pb)
    }    
    
    if(type=="Et"){
        if(any(is.na(dat))){
            warning("Warning: expectation could not be computed in all cases, see attr( . ,'empirical') to evaluate the scope of this problem",immediate.=TRUE)
        }
        predictmat <- colMeans(dat,na.rm=TRUE)
        attr(predictmat,"empirical") <- t(dat)  
    }
    
    if(length(indx==1)&(type=="hazard"|type=="survival"|type=="density")){
        predictmat <- toplot
    }
       
    
    return(list(t=t,predict=predictmat))
    
}




##' priorposterior function
##'
##' A function to produce plots of the prior (which shows as a red line) and posterior (showing as a histogram)
##'
##' @param x an object inheriting class mcmcspatsurv 
##' @param breaks see ?hist 
##' @param ylab optional y label 
##' @param main optional title 
##' @param pause logical whether to pause between plots, the default is TRUE
##' @param ... other arguments passed to the hist function
##' @return plots of the prior (red line) and posterior (histogram).
##' @seealso \link{print.mcmcspatsurv}, \link{quantile.mcmcspatsurv}, \link{summary.mcmcspatsurv}, \link{vcov.mcmcspatsurv}, 
##' \link{frailtylag1}, \link{spatialpars}, \link{hazardpars}, \link{fixedpars}, \link{randompars},
##' \link{baselinehazard}, \link{predict.mcmcspatsurv}, \link{posteriorcov}, \link{MCE},
##' \link{hazardexceedance}
##' @export

priorposterior <- function(x,breaks=30,ylab="Density",main="",pause=TRUE,...){
    nbeta <- ncol(x$betasamp)
    nomega <- ncol(x$omegasamp)
    neta <- ncol(x$etasamp)
    
    ######################
    # beta
    
    pmean <- x$priors$betaprior$mean
    psd <- x$priors$betaprior$sd
    if(length(pmean)==1){
        pmean <- rep(pmean,ncol(x$betasamp))
    }
    if(length(psd)==1){
        psd <- rep(psd,ncol(x$betasamp))
    }     
    
    for(i in 1:nbeta){
        h <- hist(x$betasamp[,i],ylab=ylab,xlab=colnames(x$betasamp)[i],breaks=breaks,freq=FALSE,main=main,...)
        xrg <- range(h$breaks)
        r <- seq(xrg[1],xrg[2],length.out=1000)               
        lines(r,dnorm(r,mean=pmean[i],sd=psd[i]),col="red",lwd=2)       
        if(pause){
            cat("[press [enter] to continue]")
            scan(n=1,quiet=TRUE)
        }
    }
    
    
    ######################
    # omega
    
    pmean <- x$priors$omegaprior$mean
    psd <- x$priors$omegaprior$sd
    if(length(pmean)==1){
        pmean <- rep(pmean,ncol(x$omegasamp))
    }
    if(length(psd)==1){
        psd <- rep(psd,ncol(x$omegasamp))
    }    
      
    samp <- x$omegasamp
    if(ncol(samp)>1){
        samp <- t(apply(samp,1,x$control$omegatrans))
    }
    else{
        samp <- t(t(apply(samp,1,x$control$omegatrans)))
    }
    colnames(samp) <- distinfo(x$dist)()$parnames
    for(i in 1:nomega){        
        h <- hist(samp[,i],ylab=ylab,xlab=paste("Transformed",colnames(x$omegasamp)[i]),breaks=breaks,freq=FALSE,main=main,...)
        xrg <- range(h$breaks)
        r <- seq(xrg[1],xrg[2],length.out=1000)
        lines(r,dnorm(r,mean=pmean[i],sd=psd[i]),col="red",lwd=2)       
        if(pause){
            cat("[press [enter] to continue]")
            scan(n=1,quiet=TRUE)
        }
    }
    
    ######################
    # eta 
    
    pmean <- x$priors$etaprior$mean
    psd <- x$priors$etaprior$sd
    if(length(pmean)==1){
        pmean <- rep(pmean,ncol(x$etasamp))
    }
    if(length(psd)==1){
        psd <- rep(psd,ncol(x$etasamp))
    }       
    
    samp <- x$etasamp
     if(ncol(samp)>1){
        samp <- t(apply(samp,1,x$cov.model$trans))
    }
    else{
        samp <- t(t(apply(samp,1,x$cov.model$trans)))
    }
    colnames(samp) <- distinfo(x$dist)()$parnames 
    for(i in 1:neta){
        h <- hist(samp[,i],ylab=ylab,xlab=paste("Transformed",colnames(x$etasamp)[i]),breaks=breaks,freq=FALSE,main=main,...)
        xrg <- range(h$breaks)
        r <- seq(xrg[1],xrg[2],length.out=1000)
        lines(r,dnorm(r,mean=pmean[i],sd=psd[i]),col="red",lwd=2)
        if(pause){
            cat("[press [enter] to continue]")
            scan(n=1,quiet=TRUE)
        }
    }
}



##' posteriorcov function
##'
##' A function to produce a plot of the posterior covariance function with upper and lower quantiles.
##'
##' @param x an object of class mcmcspatsurv 
##' @param probs vector of probabilities to be fed to quantile function  
##' @param rmax  maximum distance in space to compute this distance up to
##' @param n the number of points at which to evaluate the posterior covariance.
##' @param plot whether to plot the result
##' @param ... other arguments to be passed to matplot function 
##' @return produces a plot of the posterior spatial covariance function.
##' @seealso \link{print.mcmcspatsurv}, \link{quantile.mcmcspatsurv}, \link{summary.mcmcspatsurv}, \link{vcov.mcmcspatsurv}, 
##' \link{frailtylag1}, \link{spatialpars}, \link{hazardpars}, \link{fixedpars}, \link{randompars},
##' \link{baselinehazard}, \link{predict.mcmcspatsurv}, \link{priorposterior}, \link{MCE},
##' \link{hazardexceedance}
##' @export

posteriorcov <- function(x,probs=c(0.025,0.5,0.975),rmax=NULL,n=100,plot=TRUE,...){
    nr <- nrow(x$etasamp)
    nc <- ncol(x$etasamp)
    
    rmaxx <- 0.25*sum(apply(bbox(x$data),1,diff))/2 # approx 1/4 of mean length of observation window
    if(!is.null(rmax)){
        rmaxx <- rmax
    }
    
    r <- seq(0,rmaxx,length.out=n)
    covs <- t(apply(x$etasamp,1,function(pp){x$cov.model$eval(r,pars=pp)}))  
    qts <- t(apply(covs,2,quantile,probs=probs))
    
    rownames(qts) <- r
    
    if(plot){
        if(length(probs)==3){
            matplot(r,qts,type="l",col=c("purple","black","blue"),lty=c("dashed","solid","dashed"),xlab="Distance",ylab="Covariance")
            legend("topright",lty=c("dashed","solid","dashed"),col=rev(c("purple","black","blue")),legend=rev(probs))
        }
        else{
            matplot(r,qts,type="l",xlab="Distance",ylab="Covariance")
        }
    }
    
    return(list(r=r,qts=qts)) 
}




##' MCE function
##'
##' A function to compute Monte Carlo expectations from an object inheriting class mcmcspatsurv 
##'
##' @param object an object inheriting class mcmcspatsurv 
##' @param fun a function with arguments beta, omega, eta and Y 
##' @return the Monte Carlo mean of the function over the posterior.
##' @seealso \link{print.mcmcspatsurv}, \link{quantile.mcmcspatsurv}, \link{summary.mcmcspatsurv}, \link{vcov.mcmcspatsurv}, 
##' \link{frailtylag1}, \link{spatialpars}, \link{hazardpars}, \link{fixedpars}, \link{randompars},
##' \link{baselinehazard}, \link{predict.mcmcspatsurv}, \link{priorposterior}, \link{posteriorcov}, 
##' \link{hazardexceedance}
##' @export

MCE <- function(object,fun){
    nits <- nrow(object$betasamp)
    
    result <- lapply(1:nits,function(i){fun(beta=object$betasamp[i,],omega=object$omegasamp[i,],eta=object$etasamp[i,],Y=object$Ysamp[i,])})
    
    return((1/nits)*Reduce('+',result))
}




##' hazardexceedance function
##'
##' A function to compute exceedance probabilities for the spatially correlated frailties.
##'
##' @param threshold vector of thresholds 
##' @param direction default is "upper" which will calculate P(Y>threshold), alternative is "lower", which will calculate P(Y<threshold)
##' @return a function that can be passed to the function MCE in order to compute the exceedance probabilities
##' @seealso \link{print.mcmcspatsurv}, \link{quantile.mcmcspatsurv}, \link{summary.mcmcspatsurv}, \link{vcov.mcmcspatsurv}, 
##' \link{frailtylag1}, \link{spatialpars}, \link{hazardpars}, \link{fixedpars}, \link{randompars},
##' \link{baselinehazard}, \link{predict.mcmcspatsurv}, \link{priorposterior}, \link{posteriorcov}, \link{MCE},
##' @export

hazardexceedance <- function(threshold,direction="upper"){
    fun <- function(beta,omega,eta,Y){
        EY <- exp(Y)
        d <- length(Y)
        len <- length(threshold)
        A <- matrix(NA,len,d)
        
        for(i in 1:len){
            if(direction=="upper"){
                A[i,] <- as.numeric(EY>threshold[i])
            }
            else{
                A[i,] <- as.numeric(EY<threshold[i])
            }
        }
        return(A)        
    }
    attr(fun,"threshold") <- threshold
    attr(fun,"direction") <- direction
    return(fun)
}
