##' proposalvariance.exp function
##'
##' A function to compute an approximate proposal variance for the MALA algorithm for a parametric proportional hazards model with baseline hazard
##' derived from the exponential distribution. 
##'
##' @param X design matrix 
##' @param delta censoring indicator 
##' @param tm observed times
##' @param betahat estimate of model parameter beta 
##' @param omegahat estimate of model parameter omega
##' @param Yhat estimate of the latent field  
##' @param priors the priors, an object of class 'mcmcPriors', see ?mcmcPriors
##' @param covmodel an object of class 'covmodel', see ?covmodel
##' @param u vector of distances between points 
##' @return estimates of eta, gamma and a proposal variance matrixc for use in the MALA algorithm
##' @export

proposalvariance.exp <- function(X,delta,tm,betahat,omegahat,Yhat,priors,covmodel,u){
    
    n <- length(tm)
    lenbeta <- length(betahat)
    lenomega <- length(omegahat)
    leneta <- 2
    lenY <- length(Yhat)
    npars <- lenbeta + lenomega + leneta + lenY
    
    sigma <- matrix(0,npars,npars)
    
    # eta
    logpost <- function(eta,tm,delta,X,beta,omega,Y,priors,covmodel,u){
        sigmainv <- solve(matrix(getcov(u=u,sigma=exp(eta[1]),phi=exp(eta[2]),model=covmodel$model,pars=covmodel$pars),n,n))
        cholsigmainv <- t(chol(sigmainv))
        gamma <- cholsigmainv%*%(Y+exp(eta[1])^2/2)                    
        n <- nrow(X)
        Xbeta <- X%*%beta        
        priorcontrib <- -(1/2)*sum(gamma^2) + do.call(priors$call,args=list(beta=beta,omega=omega,eta=eta,priors=priors))
        stuff <- Xbeta + Y + omega
        expstuff <- exp(stuff)
        logpost <- sum(log(diag(cholsigmainv))) + sum(delta*(stuff)-expstuff*tm) + priorcontrib # first term, sum(log(diag(cholsigmainv))), is the Jacobian
        return(logpost)
    }
    ngrid <- 20
    if(length(priors$etaprior$mean)==1){
        xseq <- seq(priors$etaprior$mean-1.96*priors$etaprior$sd,priors$etaprior$mean+1.96*priors$etaprior$sd,length.out=ngrid)
        yseq <- xseq
    }
    else{
        xseq <- seq(priors$etaprior$mean[1]-1.96*priors$etaprior$sd[1],priors$etaprior$mean[1]+1.96*priors$etaprior$sd[1],length.out=ngrid)
        yseq <- seq(priors$etaprior$mean[2]-1.96*priors$etaprior$sd[2],priors$etaprior$mean[2]+1.96*priors$etaprior$sd[2],length.out=ngrid)
    }
    qa <- quadapprox(logpost,xseq=xseq,yseq=yseq,tm=tm,delta=delta,X=X,beta=betahat,omega=omegahat,Y=Yhat,priors=priors,covmodel=covmodel,u=u)
    matr <- 0.4*qa$curvature
    etahat <- qa$max
    
    # entry for eta in propossal covariance
    sigma[(lenbeta+lenomega+1):(lenbeta+lenomega+leneta),(lenbeta+lenomega+1):(lenbeta+lenomega+leneta)] <- matr    
    
    #estimate of gamma
    Sigma <- matrix(getcov(u=u,sigma=exp(etahat[1]),phi=exp(etahat[2]),model=covmodel$model,pars=covmodel$pars),n,n)
    covinv <- solve(Sigma)
    cholcovinv <- t(chol(covinv))
    gammahat <- cholcovinv%*%(Yhat+exp(etahat[1])^2/2)  
    
    cholSigma <- t(chol(Sigma))    
    #mu <- -sd(gammahat)^2/2
    #Y <- mu + cholSigma%*%gammahat
    
    deriv <- do.call(priors$derivative,args=list(beta=betahat,omega=omegahat,eta=etahat,priors=priors))

    N <- nrow(X)
    d <- ncol(X)
    
    Xbetahat <- X%*%betahat
    thing <- exp(Xbetahat + Yhat + omegahat)
    
    # beta and omega
    for(k in 1:lenbeta){
        for(l in 1:lenbeta){
            sigma[k,l] <- sum(-X[,k]*X[,l]*thing*tm) + as.numeric(k==l)*deriv$deriv2[k]
        }
        sigma[k,lenbeta+1] <- sigma[d+1,k] <- sum(-X[,k]*thing*tm) 
    }
    sigma[lenbeta+1,lenbeta+1] <- -sum(thing*tm) + deriv$deriv2[lenbeta+1] 
    
    # gamma
    diag(sigma)[(lenbeta+lenomega+leneta+1):npars] <- -sum(diag(cholSigma^2)%*%(thing*tm)) - 1  # -1 comes from prior
    
    return(list(etahat=etahat,gammahat=gammahat,sigma=solve(-sigma))) 
}



##' estimateY.exp function
##'
##' A function to estimate Y assuming a parametric proportional hazards model with baseline hazard derived from the exponential distribution.
##'
##' @param X the design matrix
##' @param betahat an estimate of beta
##' @param omegahat an estimate of omega
##' @param tm vector of observed times
##' @param delta censoring indicator
##' @param u vecor of distances bettween points
##' @param covmodel an object of class 'covmodel', see ?covmodel 
##' @return an estimate of the latent Gaussian field
##' @export

estimateY.exp <- function(X,betahat,omegahat,tm,delta,u,covmodel){
    
    tsubs <- tm
    for(i in 1:length(tsubs)){
        if(delta[i]==0){
            tpot <- tsubs[tsubs>tsubs[i]] # potential t
            if(length(tpot)==0){
                next # leave tsubs[i] alone 
            }
            else{
                tsubs[i] <- sample(tpot,1) # sample from empirical distribution (ignoring covariates)
            }
        }
    }
    Y <- -log(tsubs) - X%*%betahat - omegahat # greedy estimate of Y (maximise individual contributions to log-likelihood) ... note log(delta) is now omitted

    return(Y)    
}
