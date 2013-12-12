##' transformestimates.exp function
##'
##' A function to transform estimates of the parameters of the exponential baseline hazard function, so they are commensurate with R's inbuilt density functions.
##'
##' @param x a vector of paramters 
##' @return the transformed parameters. For the exponential model this is just the identity.
##' @export
transformestimates.exp <- function(x){
    return(x) # note this is logged later for use in the MCMC
}



##' invtransformestimates.exp function
##'
##' A function to transform estimates of the parameters of the exponential baseline hazard function, so they are commensurate with R's inbuilt density functions.
##'
##' @param x a vector of paramters 
##' @return the transformed parameters. For the exponential model this is just the identity.
##' @export
invtransformestimates.exp <- function(x){
    return(x)
}




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
##' @param control list of control parameters, see ?inference.control
##' @return estimates of eta, gamma and a proposal variance matrixc for use in the MALA algorithm
##' @export
proposalvariance.exp <- function(X,delta,tm,betahat,omegahat,Yhat,priors,covmodel,u,control){
     
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
        gamma <- cholsigmainv%*%(Y-exp(eta[1])^2/2)                    
        n <- nrow(X)
        Xbeta <- X%*%beta        
        priorcontrib <- -(1/2)*sum(gamma^2) + do.call(priors$call,args=list(beta=beta,omega=omega,eta=eta,priors=priors))
        stuff <- Xbeta + Y + omega
        expstuff <- exp(stuff)
        logpost <- sum(delta*(stuff)-expstuff*tm) + priorcontrib # first term, sum(log(diag(cholsigmainv))), is the Jacobian
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
    matr <- qa$curvature
    etahat <- qa$max
    
    # entry for eta in propossal covariance
    sigma[(lenbeta+lenomega+1):(lenbeta+lenomega+leneta),(lenbeta+lenomega+1):(lenbeta+lenomega+leneta)] <- matr    
    
    #estimate of gamma
    Sigma <- matrix(getcov(u=u,sigma=exp(etahat[1]),phi=exp(etahat[2]),model=covmodel$model,pars=covmodel$pars),n,n)
    covinv <- solve(Sigma)
    cholcovinv <- t(chol(covinv))
    gammahat <- cholcovinv%*%(Yhat-exp(etahat[1])^2/2)  
    
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
    diag(sigma)[(lenbeta+lenomega+leneta+1):npars] <- -sum((diag(cholSigma)^2)%*%(thing*tm)) - 1  # -1 comes from prior  
    
    return(list(etahat=etahat,gammahat=gammahat,sigma=solve(-sigma))) 
}


##' proposalvariance.exp.gridded function
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
##' @param control list of control parameters, see ?inference.control
##' @return estimates of eta, gamma and a proposal variance matrixc for use in the MALA algorithm
##' @export
proposalvariance.exp.gridded <- function(X,delta,tm,betahat,omegahat,Yhat,priors,covmodel,u,control){
      
    Ygrid <- gridY(Y=Yhat,control=control)
    
    
    n <- length(tm)
    lenbeta <- length(betahat)
    lenomega <- length(omegahat)
    leneta <- 2
    lenY <- length(Ygrid)
    npars <- lenbeta + lenomega + leneta + lenY
    
    #sigma <- matrix(0,npars,npars)
    sigma <- Matrix(0,npars,npars)

    
    # eta
    logpost <- function(eta,tm,delta,X,beta,omega,Ygrid,priors,covmodel,u,control){
        
        covbase <- matrix(getcov(u=u,sigma=exp(eta[1]),phi=exp(eta[2]),model=covmodel$model,pars=covmodel$pars),control$Mext,control$Next)
        rootQeigs <- sqrt(1/Re(fft(covbase)))
        
        gamma <- GammafromY(Ygrid,rootQeigs=rootQeigs,mu=-(exp(eta[1]))^2/2)   
                          
        n <- nrow(X)
        Xbeta <- X%*%beta        
        priorcontrib <- -(1/2)*sum(gamma^2) + do.call(priors$call,args=list(beta=beta,omega=omega,eta=eta,priors=priors))
        stuff <- Xbeta + Ygrid[control$idx] + omega
        expstuff <- exp(stuff)
        logpost <- sum(delta*(stuff)-expstuff*tm) + priorcontrib
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
    qa <- quadapprox(logpost,xseq=xseq,yseq=yseq,tm=tm,delta=delta,X=X,beta=betahat,omega=omegahat,Ygrid=Ygrid,priors=priors,covmodel=covmodel,u=u,control=control)
    matr <- qa$curvature
    etahat <- qa$max
    
    # entry for eta in propossal covariance
    sigma[(lenbeta+lenomega+1):(lenbeta+lenomega+leneta),(lenbeta+lenomega+1):(lenbeta+lenomega+leneta)] <- matr    
    
    #estimate of gamma
    covbase <- matrix(getcov(u=u,sigma=exp(etahat[1]),phi=exp(etahat[2]),model=covmodel$model,pars=covmodel$pars),control$Mext,control$Next)
    rootQeigs <- sqrt(1/Re(fft(covbase)))
    invrootQeigs <- 1/rootQeigs
    gammahat <- GammafromY(Ygrid,rootQeigs=rootQeigs,mu=-(exp(etahat[1]))^2/2)     
        
    deriv <- do.call(priors$derivative,args=list(beta=betahat,omega=omegahat,eta=etahat,priors=priors))

    N <- nrow(X)
    d <- ncol(X)
    
    Xbetahat <- X%*%betahat
    thing <- exp(Xbetahat + Ygrid[control$idx] + omegahat)
    
    # beta and omega
    for(k in 1:lenbeta){
        for(l in 1:lenbeta){
            sigma[k,l] <- sum(-X[,k]*X[,l]*thing*tm) + as.numeric(k==l)*deriv$deriv2[k]
        }
        sigma[k,lenbeta+1] <- sigma[d+1,k] <- sum(-X[,k]*thing*tm) 
    }
    sigma[lenbeta+1,lenbeta+1] <- -sum(thing*tm) + deriv$deriv2[lenbeta+1] 
    
    # gamma
    ident <- matrix(0,control$Mext,control$Next)
    ident[1,1] <- 1 # gives the identity matrix in block circulant form
    cholmat <- Re((1/(control$Mext*control$Next))*fft(invrootQeigs*fft(ident,inverse=TRUE)))
    cholmat2 <- cholmat^2
    
    thingaug <- matrix(0,control$Mext,control$Next)
    thingtemp <- sapply(control$uqidx,function(i){sum(thing[control$idx==i])})
    thingaug[control$uqidx] <- thingtemp

    matidx <- (lenbeta+lenomega+leneta+1):npars
    matidx <- matrix(matidx,nrow=length(matidx),ncol=2) 
    
    jobby <- sapply(control$uqidx,function(i){return(sum(thing[control$idx==i]*tm[control$idx==i]))})
    sigma[matidx] <- -1 # -1 comes from prior 
    matidx1 <- lenbeta+lenomega+leneta+control$uqidx
    matidx1 <- matrix(matidx1,nrow=length(matidx1),ncol=2) 
    sigma[matidx1] <- sigma[matidx1] - length(jobby)*cholmat2[1,1]*jobby    
    
    sigma <- (-1) * sigma # variance is inverse of observed information    
    
    sigmaret <- Matrix(0,npars,npars)
    sigmaret[1:(lenbeta+lenomega+leneta),1:(lenbeta+lenomega+leneta)] <- solve(as.matrix(sigma[1:(lenbeta+lenomega+leneta),1:(lenbeta+lenomega+leneta)]))
    sigmaret[matidx] <- 1/sigma[matidx]   
    
    return(list(etahat=etahat,gammahat=gammahat,sigma=sigmaret)) 
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
##' @return an estimate of the latent Gaussian field
##' @export
estimateY.exp <- function(X,betahat,omegahat,tm,delta){
    
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
