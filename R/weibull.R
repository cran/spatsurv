##' transformestimates.weibull function
##'
##' A function to transform estimates of the parameters of the weibull baseline hazard function, so they are commensurate with R's inbuilt density functions.
##'
##' @param x a vector of paramters
##' @return the transformed parameters. For the weibull model, this transforms 'shape' 'scale' (see ?dweibull) to 'alpha' and 'lambda' for the MCMC
##' @export
transformestimates.weibull <- function(x){
    a <- x[1] # shape
    b <- x[2] # scale
    alpha <- a
    lambda <- (1/b)^a

    ans <- c(alpha=alpha,lambda=lambda) # note this is logged later for use in the MCMC
    
    return(ans)
}



##' invtransformestimates.weibull function
##'
##' A function to transform estimates of the parameters of the weibull baseline hazard function, so they are commensurate with R's inbuilt density functions.
##'
##' @param x a vector of paramters
##' @return the transformed parameters. For the weibull model, this is the back-transform from 'alpha' and 'lambda' to 'shape' 'scale' (see ?dweibull).
##' @export
invtransformestimates.weibull <- function(x){

    alpha <- x[1]
    lambda <- x[2]
    
    shape <- alpha
    scale <- exp((-1/alpha)*log(lambda))

    ans <- c(shape=shape,scale=scale)    
    
    return(ans)
}

##' proposalvariance.weibull function
##'
##' A function to compute an approximate proposal variance for the MALA algorithm for a parametric proportional hazards model with baseline hazard
##' derived from the Weibull distribution.
##'
##' @param X design matrix 
##' @param delta censoring indicator 
##' @param tm observed times
##' @param betahat estimate of model parameter beta 
##' @param omegahat estimate of model parameter omega
##' @param Yhat estimate of the latent field  
##' @param priors the priors, an object of class 'mcmcPriors', see ?mcmcPriors
##' @param cov.model an object of class 'covmodel', see ?covmodel
##' @param u vector of distances between points
##' @param control list of control parameters, see ?inference.control
##' @return estimates of eta, gamma and a proposal variance matrixc for use in the MALA algorithm
##' @export
proposalvariance.weibull <- function(X,delta,tm,betahat,omegahat,Yhat,priors,cov.model,u,control){
    
    n <- length(tm)
    lenbeta <- length(betahat)
    lenomega <- length(omegahat)
    leneta <- getleneta(cov.model)
    lenY <- length(Yhat)
    npars <- lenbeta + lenomega + leneta + lenY
    
    sigma <- matrix(0,npars,npars)
    
    # eta
    logpost <- function(eta,tm,delta,X,beta,omega,Y,priors,cov.model,u){
        #sigmainv <- solve(matrix(getcov(u=u,sigma=exp(eta[1]),phi=exp(eta[2]),model=cov.model$model,pars=cov.model$pars),n,n))
        sigmainv <- solve(matrix(EvalCov(cov.model=cov.model,u=u,parameters=eta),n,n))
        cholsigmainv <- t(chol(sigmainv))
        gamma <- cholsigmainv%*%(Y-exp(eta[1])^2/2)                    
        
        alpha <- exp(omega[1])
        lambda <- exp(omega[2])
    
        n <- nrow(X)
        Xbeta <- X%*%beta
        #sigma <- matrix(getcov(u=u,sigma=exp(eta[1]),phi=exp(eta[2]),model=cov.model$model,pars=cov.model$pars),n,n)
        #cholsigma <- t(chol(sigma))
        priorcontrib <- -(1/2)*sum(gamma^2) + do.call(priors$call,args=list(beta=beta,omega=omega,eta=eta,priors=priors))
        #Y <- -eta[1]^2/2 + cholsigma%*%gamma
        stuff <- Xbeta + Y
        expstuff <- exp(stuff)
    
        logpost <- sum(delta*(stuff + log(lambda) + log(alpha) + (alpha-1)*log(tm))-expstuff*lambda*tm^alpha) + priorcontrib # first term, sum(log(diag(cholsigmainv))), is the Jacobian
    
        return(logpost)
    }
    
    npts <- 20
    if(leneta>=3){
        npts <- 10
    }
    rgs <- getparranges(priors=priors,leneta=leneta)   
    qa <- QuadApprox(logpost,npts=npts,argRanges=rgs,tm=tm,delta=delta,X=X,beta=betahat,omega=omegahat,Y=Yhat,priors=priors,cov.model=cov.model,u=u)    
    
    matr <- qa$curvature
    etahat <- qa$max
    
    # entry for eta in proposal covariance
    sigma[(lenbeta+lenomega+1):(lenbeta+lenomega+leneta),(lenbeta+lenomega+1):(lenbeta+lenomega+leneta)] <- matr    
    
    #estimate of gamma
    #Sigma <- matrix(getcov(u=u,sigma=exp(etahat[1]),phi=exp(etahat[2]),model=cov.model$model,pars=cov.model$pars),n,n)
    Sigma <- matrix(EvalCov(cov.model=cov.model,u=u,parameters=etahat),n,n)
    covinv <- solve(Sigma)
    cholcovinv <- t(chol(covinv))
    gammahat <- cholcovinv%*%(Yhat-exp(etahat[1])^2/2)  
    
    cholSigma <- t(chol(Sigma))    
    #mu <- -sd(gammahat)^2/2
    #Y <- mu + cholSigma%*%gammahat
    
    deriv <- do.call(priors$derivative,args=list(beta=betahat,omega=omegahat,eta=etahat,priors=priors))

    N <- nrow(X)
    d <- ncol(X)
    
    alphahat <- exp(omegahat[1])
    lambdahat <- exp(omegahat[2])
    
    alphaJacobian <- alphahat
    alpha2Jacobian <- alphahat
    
    lambdaJacobian <- lambdahat
    lambda2Jacobian <- lambdahat    
    
    Xbetahat <- X%*%betahat
    thing <- exp(Xbetahat + Yhat)
    
    # beta and omega
    for(k in 1:lenbeta){
        for(l in 1:lenbeta){
            sigma[k,l] <- sum(-X[,k]*X[,l]*thing*lambdahat*tm^alphahat) + as.numeric(k==l)*deriv$deriv2[k]
        }
        sigma[k,d+1] <- sigma[d+1,k] <- alphaJacobian * sum(-X[,k]*thing*lambdahat*log(tm)*tm^alphahat)
        sigma[k,d+2] <- sigma[d+2,k] <- lambdaJacobian * sum(-X[,k]*thing*tm^alphahat) 
    }
    sigma[lenbeta+1,lenbeta+2] <- sigma[lenbeta+2,lenbeta+1] <- alphaJacobian * lambdaJacobian * sum(-thing*log(tm)*tm^alphahat)
    sigma[lenbeta+1,lenbeta+1] <- alphaJacobian^2*sum(-delta/(alphahat^2)-thing*lambdahat*(log(tm))^2*tm^alphahat) + alpha2Jacobian*sum(delta*(1/alphahat+log(tm))-thing*lambdahat*log(tm)*tm^alphahat)
    sigma[lenbeta+2,lenbeta+2] <- lambdaJacobian^2 * sum(-delta/lambdahat^2) + lambda2Jacobian*sum(delta/lambdahat-thing*tm^alphahat)  
    
    # gamma
    diag(sigma)[(lenbeta+lenomega+leneta+1):npars] <- -sum((diag(cholSigma)^2)*(thing*lambdahat*tm^alphahat)) - 1 # -1 comes from prior 
    
    return(list(etahat=etahat,gammahat=gammahat,sigma=solve(-sigma))) 
}



##' proposalvariance.weibull.gridded function
##'
##' A function to compute an approximate proposal variance for the MALA algorithm for a parametric proportional hazards model with baseline hazard
##' derived from the Weibull distribution.
##'
##' @param X design matrix 
##' @param delta censoring indicator 
##' @param tm observed times
##' @param betahat estimate of model parameter beta 
##' @param omegahat estimate of model parameter omega
##' @param Yhat estimate of the latent field  
##' @param priors the priors, an object of class 'mcmcPriors', see ?mcmcPriors
##' @param cov.model an object of class 'covmodel', see ?covmodel
##' @param u vector of distances between points
##' @param control list of control parameters, see ?inference.control
##' @return estimates of eta, gamma and a proposal variance matrixc for use in the MALA algorithm
##' @export
proposalvariance.weibull.gridded <- function(X,delta,tm,betahat,omegahat,Yhat,priors,cov.model,u,control){
    
    Ygrid <- gridY(Y=Yhat,control=control)    
    
    n <- length(tm)
    lenbeta <- length(betahat)
    lenomega <- length(omegahat)
    leneta <- getleneta(cov.model)
    lenY <- length(Ygrid)
    npars <- lenbeta + lenomega + leneta + lenY
    
    #sigma <- matrix(0,npars,npars)
    sigma <- Matrix(0,npars,npars)
    
    # eta
    logpost <- function(eta,tm,delta,X,beta,omega,Ygrid,priors,cov.model,u,control){

        #covbase <- matrix(getcov(u=u,sigma=exp(eta[1]),phi=exp(eta[2]),model=cov.model$model,pars=cov.model$pars),control$Mext,control$Next)
        covbase <- matrix(EvalCov(cov.model=cov.model,u=u,parameters=eta),control$Mext,control$Next)
        rootQeigs <- sqrt(1/Re(fft(covbase)))
        
        #browser()
        pars <- sapply(1:length(eta),function(i){cov.model$itrans[[i]](eta[i])})
        ymean <- -pars[which(cov.model$parnames=="sigma")]^2/2
        gamma <- GammafromY(Ygrid,rootQeigs=rootQeigs,mu=ymean)
        #gamma <- GammafromY(Ygrid,rootQeigs=rootQeigs,mu=-(exp(eta[1]))^2/2)                      
        
        alpha <- exp(omega[1])
        lambda <- exp(omega[2])
    
        n <- nrow(X)
        Xbeta <- X%*%beta
        priorcontrib <- -(1/2)*sum(gamma^2) + do.call(priors$call,args=list(beta=beta,omega=omega,eta=eta,priors=priors))
        stuff <- Xbeta + Ygrid[control$idx]
        expstuff <- exp(stuff)
    
        logpost <- sum(delta*(stuff + log(lambda) + log(alpha) + (alpha-1)*log(tm))-expstuff*lambda*tm^alpha) + priorcontrib # first term, sum(log(diag(cholsigmainv))), is the Jacobian
    
        return(logpost)
    }
    
    npts <- 20
    if(leneta>=3){
        npts <- 10
    }
    rgs <- getparranges(priors=priors,leneta=leneta)   
    qa <- QuadApprox(logpost,npts=npts,argRanges=rgs,tm=tm,delta=delta,X=X,beta=betahat,omega=omegahat,Ygrid=Ygrid,priors=priors,cov.model=cov.model,u=u,control=control)    
    
    matr <- qa$curvature
    etahat <- qa$max
    
    # entry for eta in proposal covariance
    sigma[(lenbeta+lenomega+1):(lenbeta+lenomega+leneta),(lenbeta+lenomega+1):(lenbeta+lenomega+leneta)] <- matr    
    
    #estimate of gamma
    covbase <- matrix(EvalCov(cov.model=cov.model,u=u,parameters=etahat),control$Mext,control$Next)
    rootQeigs <- sqrt(1/Re(fft(covbase)))
    invrootQeigs <- 1/rootQeigs
    #gammahat <- GammafromY(Ygrid,rootQeigs=rootQeigs,mu=-(exp(etahat[1]))^2/2)    
    
    deriv <- do.call(priors$derivative,args=list(beta=betahat,omega=omegahat,eta=etahat,priors=priors))

    N <- nrow(X)
    d <- ncol(X)
    
    alphahat <- exp(omegahat[1])
    lambdahat <- exp(omegahat[2])
    
    alphaJacobian <- alphahat
    alpha2Jacobian <- alphahat
    
    lambdaJacobian <- lambdahat
    lambda2Jacobian <- lambdahat    
    
    Xbetahat <- X%*%betahat
    thing <- exp(Xbetahat + Ygrid[control$idx])
    
    # beta and omega
    for(k in 1:lenbeta){
        for(l in 1:lenbeta){
            sigma[k,l] <- sum(-X[,k]*X[,l]*thing*lambdahat*tm^alphahat) + as.numeric(k==l)*deriv$deriv2[k]
        }
        sigma[k,d+1] <- sigma[d+1,k] <- alphaJacobian * sum(-X[,k]*thing*lambdahat*log(tm)*tm^alphahat)
        sigma[k,d+2] <- sigma[d+2,k] <- lambdaJacobian * sum(-X[,k]*thing*tm^alphahat) 
    }
    sigma[lenbeta+1,lenbeta+2] <- sigma[lenbeta+2,lenbeta+1] <- alphaJacobian * lambdaJacobian * sum(-thing*log(tm)*tm^alphahat)
    sigma[lenbeta+1,lenbeta+1] <- alphaJacobian^2*sum(-delta/(alphahat^2)-thing*lambdahat*(log(tm))^2*tm^alphahat) + alpha2Jacobian*sum(delta*(1/alphahat+log(tm))-thing*lambdahat*log(tm)*tm^alphahat)
    sigma[lenbeta+2,lenbeta+2] <- lambdaJacobian^2 * sum(-delta/lambdahat^2) + lambda2Jacobian*sum(delta/lambdahat-thing*tm^alphahat)  
    
    # gamma
    #diag(sigma)[(lenbeta+lenomega+leneta+1):npars] <- -sum((diag(cholSigma)^2)*(thing*lambdahat*tm^alphahat)) - 1 # -1 comes from prior 
    
    #browser()    
    
    ident <- matrix(0,control$Mext,control$Next)
    ident[1,1] <- 1 # gives the identity matrix in block circulant form
    cholmat <- Re((1/(control$Mext*control$Next))*fft(invrootQeigs*fft(ident,inverse=TRUE)))
    cholmat2 <- cholmat^2
    
    thingaug <- matrix(0,control$Mext,control$Next)
    thingtemp <- sapply(control$uqidx,function(i){sum(thing[control$idx==i])})
    thingaug[control$uqidx] <- thingtemp

    matidx <- (lenbeta+lenomega+leneta+1):npars
    matidx <- matrix(matidx,nrow=length(matidx),ncol=2) 
    
    jobby <- sapply(control$uqidx,function(i){return(sum(thing[control$idx==i]*lambdahat*tm[control$idx==i]^alphahat))})
    sigma[matidx] <- -1 # -1 comes from prior 
    matidx1 <- lenbeta+lenomega+leneta+control$uqidx
    matidx1 <- matrix(matidx1,nrow=length(matidx1),ncol=2) 
    sigma[matidx1] <- sigma[matidx1] - length(jobby)*cholmat2[1,1]*jobby    
    
    sigma <- (-1) * sigma # variance is inverse of observed information    
    
    sigmaret <- Matrix(0,npars,npars)
    sigmaret[1:(lenbeta+lenomega+leneta),1:(lenbeta+lenomega+leneta)] <- solve(as.matrix(sigma[1:(lenbeta+lenomega+leneta),1:(lenbeta+lenomega+leneta)]))
    sigmaret[matidx] <- 1/sigma[matidx] 
    
    return(list(etahat=etahat,sigma=sigmaret)) 
}



##' estimateY.weibull function
##'
##' A function to estimate Y assuming a parametric proportional hazards model with baseline hazard derived from the Weibull distribution.
##'
##' @param X the design matrix
##' @param betahat an estimate of beta
##' @param omegahat an estimate of omega
##' @param tm vector of observed times
##' @param delta censoring indicator
##' @return an estimate of the latent Gaussian field
##' @export
estimateY.weibull <- function(X,betahat,omegahat,tm,delta){

    alpha <- exp(omegahat[1])
    lambda <- exp(omegahat[2])
    
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
    
    Y <- -X%*%betahat - log(lambda) - alpha*log(tsubs) # greedy estimate of Y (maximise individual contributions to log-likelihood) ... note log(delta) is now omitted  

    return(Y)    
}
