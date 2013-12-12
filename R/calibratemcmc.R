##' proposalvariance_nospat function
##'
##' A function to compute an approximate scaling matrix for use with MCMC algorithms. Works for models with baseline hazard derived from the
##' exponential survival model and assumes no frailties i.e. this is used for a non-spatial survival analysis.
##'
##' @param X matrix of covariates 
##' @param delta censoring indicator, a vector 
##' @param tm observed event times, a vector
##' @param betahat an estimate of beta 
##' @param omegahat an estimate of omega
##' @param betapriorsd standard deviation of the prior for beta 
##' @param omegapriorsd standard deviation of the prior for omega
##' @return ...
##' @export

proposalvariance_exponential_nospat <- function(X,delta,tm,betahat,omegahat,betapriorsd,omegapriorsd){
    N <- nrow(X)
    d <- ncol(X)
    sigma <- matrix(NA,d+1,d+1)
    Xbetahat <- X%*%betahat
    expXbetahat <- exp(Xbetahat)
    
    for(k in 1:d){
        for(l in 1:d){
            sigma[k,l] <- sum(-X[,k]*X[,l]*expXbetahat*exp(omegahat)*tm) - as.numeric(k==l)/betapriorsd^2
        }
        sigma[k,d+1] <- sigma[d+1,k] <- sum(-X[,k]*expXbetahat*exp(omegahat)*tm) 
    }
    sigma[d+1,d+1] <- -sum(expXbetahat*exp(omegahat)*tm) - 1/omegapriorsd^2
    
    return(solve(-sigma))      
}




##' derivlogindepGaussianprior function
##'
##' A function to compute the first and second derivatives of the log-density assuming independent Gaussian priors for each of the parameters.
##'
##' @param beta a vector, the parameter beta
##' @param omega a vector, the parameter omega
##' @param eta a vector, the parameter eta 
##' @param priors an object of class 'mcmcPrior', see ?mcmcPrior
##' @return ...
##' @export

derivlogindepGaussianprior <- function(beta=NULL,omega=NULL,eta=NULL,priors){
    deriv1 <- c((-1/priors$betaprior$sd^2)*(beta-priors$betaprior$mean),(-1/priors$omegaprior$sd^2)*(omega-priors$omegaprior$mean),(-1/priors$etaprior$sd^2)*(eta-priors$etaprior$mean))
    sdbeta <- priors$betaprior$sd
    sdomega <- priors$omegaprior$sd
    sdeta <- priors$etaprior$sd
    if (length(priors$betaprior$sd)<length(beta)){
        sdbeta <- rep(priors$betaprior$sd,length(beta))
    }
    if (length(priors$omegaprior$sd)<length(omega)){
        sdomega <- rep(priors$omegaprior$sd,length(omega))
    }
    if (length(priors$etaprior$sd)<length(eta)){
        sdeta <- rep(priors$etaprior$sd,length(eta))
    }
    deriv2 <- c(-1/sdbeta^2,-1/sdomega^2,-1/sdeta^2)
    return(list(deriv1=deriv1,deriv2=deriv2))
}




##' QuadApprox function
##'
##' A function to compute the second derivative of a function (of several real variables) using a quadratic approximation  on a 
##' grid of points defined by the list argRanges. Also returns the local maximum. 
##'
##' @param fun a function
##' @param npts integer number of points in each direction
##' @param argRanges a list of ranges on which to construct the grid for each parameter 
##' @param plot whether to plot the quadratic approximation of the posterior (for two-dimensional parameters only)
##' @param ... other arguments to be passed to fun
##' @return a 2 by 2 matrix containing the curvature at the maximum and the (x,y) value at which the maximum occurs 
##' @export


QuadApprox <- function(fun,npts,argRanges,plot=TRUE,...){

    npar <- length(argRanges)
    vals <- lapply(argRanges,function(x){seq(x[1],x[2],length.out=npts)})
    parn <- paste("x",1:npar,sep="")
    parnames <- parn # later this will be used in the call to lm
    paridx <- as.list(rep(NA,1+2*npar+choose(npar,2))) # this will be used later to find the function maximum via a set of simultaneous equations (intercept, single, squared and mixed terms)
    paridx[(1+1:npar)] <- 1:npar
    gr <- expand.grid(vals)
    gr2 <- gr^2
    parnames <- c(parnames,paste("x",1:npar,".2",sep=""))
    paridx[(npar+2):(2*npar+1)] <- 1:npar
    grcross <- matrix(NA,nrow(gr),choose(npar,2))
    ct <- 1
    for(i in 1:(npar-1)){
        for(j in (i+1):npar){    
            grcross[,ct] <- gr[,i]*gr[,j]
            parnames <- c(parnames,paste(parn[i],parn[j],collapse="",sep=""))
            paridx[[2*npar+1+ct]] <- c(i,j)
            ct <- ct + 1
        }
    }
    partype <- c("intercept",rep("single",npar),rep("squared",npar),rep("mixed",choose(npar,2)))
    dataf <- cbind(gr,gr2,grcross)
    names(dataf) <- parnames
    
    cat("Constructing quadratic approximation to posterior (this can take some time) ...\n")
    dataf$funvals <- apply(gr,1,function(params){fun(params,...)})
    cat("Done.\n")
    
    
    if(plot){
        if(npar==2){
            image.plot(vals[[1]],vals[[2]],matrix(dataf$funvals,npts,npts),main="Function")
        }  
    }
    
    form <- paste("funvals ~",paste(parnames,collapse=" + "))
     
    mod <- lm(form,data=dataf)
    co <- coefficients(mod)
    
    if(plot){
        if(npar==2){
            image.plot(vals[[1]],vals[[2]],matrix(fitted(mod),npts,npts),main="Quadratic Approximation")
        }  
    }    
    
    # now construct matrix of second derivatives
    sigmainv <- matrix(NA,npar,npar)
    diag(sigmainv) <- 2 * co[which(partype=="squared")] # first the diagonal elements
    idx <- which(partype=="mixed") # now the off diagonals
    ct <- 1
    for(i in 1:(npar-1)){
        for(j in (i+1):npar){    
            sigmainv[i,j] <- co[idx[ct]]
            sigmainv[j,i] <- co[idx[ct]]
            ct <- ct + 1
        }
    }
    
    # lastly, create a system of simultaneous equations, Ax = b, which when solved gives the maximum
    b <- (-1) * matrix(co[which(partype=="single")],npar,1)
    A <- matrix(NA,npar,npar)
    diag(A) <- 2 * co[which(partype=="squared")]
    for(i in 1:(npar-1)){
        for(j in (i+1):npar){
            tst <- sapply(paridx,function(x){any(x==i)&any(x==j)})
            idx <- which(tst)   
            A[i,j] <- co[idx]
            A[j,i] <- co[idx]
        }
    }

    etaest <- as.vector(solve(A)%*%b) # now solve the system of simultaneous equations to get an initial guess for eta 
    
    sigmainv <- fixmatrix(sigmainv)

    return(list(max=etaest,curvature=sigmainv,mod=mod)) 
}




##' fixmatrix function
##'
##' A function to 
##'
##' @param mat X 
##' @return ...
##' @export

fixmatrix <- function(mat){
    
    mat <- (-1)*mat # since mat is curvature, the negative *should* have positive eigenvalues
    ev <- eigen(mat)$values
    if(all(ev>0)){
        return((-1)*mat)
    }
    else if(all(ev<0)){
        stop("Estimated covariance matrix for eta has all negative eigenvalues")
    }
    else{
        warning("Something is wrong with the estimated covariance matrix, fixing this using a totally ad-hoc method. This will not affect ergodicity, merely the efficiency of the chain.",immediate.=TRUE)
        cat("Fixing non positive definite covariance matrix for eta ...\n") 
    
        diag(mat) <- abs(diag(mat)) # hmmmm ....        
               
        if(all(dim(mat)==2)){
            fun <- function(x){
                tmp <- mat
                tmp[1,2] <- tmp[2,1] <- mat[1,2] / x
                posev <- abs(ev)
                ev1 <- eigen(tmp)$values
                if(!all(ev1>0)){
                    return(.Machine$double.xmax)
                }
                else{
                    df1 <- (posev[1]-ev1[1])/posev[1]
                    df2 <- (posev[2]-ev1[2])/posev[2]
                    return(df1^2+df2^2)
                }                
            }
            op <- suppressWarnings(try(optimise(fun,interval=c(0,10))))
            if(inherits(op,"try-error")){
                stop("Failed to fix negative definite matix")
            }
            ans <- mat
            ans[1,2] <- ans[2,1] <- mat[1,2] / op$minimum
                       
        }
        else{
            #browser()
            fun1 <- function(pars){
                tmp <- mat
                tmp[lower.tri(tmp)] <- tmp[lower.tri(tmp)] / pars
                tmp[upper.tri(tmp)] <- tmp[upper.tri(tmp)] / pars
                posev <- abs(ev)
                ev1 <- eigen(tmp)$values
                if(!all(ev1>0)){
                    return(.Machine$double.xmax)
                }
                else{
                    dff <- sum(((posev-ev1)/posev)^2)
                    return(dff)
                }                
            }
            op <- suppressWarnings(try(optim(par=rep(1,ncol(mat)),fn=fun1)))
            if(inherits(op,"try-error")){
                stop("Failed to fix negative definite matix")
            }
            
            ans <- mat
            ans[lower.tri(ans)] <- ans[lower.tri(ans)] / op$par
            ans[upper.tri(ans)] <- ans[upper.tri(ans)] / op$par
        }
        

        ct <- nrow(ans)        
        if(!all(eigen(ans)$values>0)){
            while(!all(eigen(ans)$values>0) & ct>0){
                ans[ct,1:(ct-1)] <- 0
                ans[1:(ct-1),ct] <- 0
                ct <- ct - 1
            }    
        }
        
        if(!all(eigen(ans)$values>0)){
            stop("Failed to fix negative definite matix")
        }
        
        ans <- (-1)*ans              
        
        return(ans)    
    } 
}



##' proposalVariance function
##'
##' A function to 
##'
##' @param X X 
##' @param surv X 
##' @param betahat X 
##' @param omegahat X 
##' @param Yhat X 
##' @param priors X 
##' @param cov.model X 
##' @param u X 
##' @param control X 
##' @return ...
##' @export

proposalVariance <- function(X,surv,betahat,omegahat,Yhat,priors,cov.model,u,control){
     
    n <- nrow(X)
    lenbeta <- length(betahat)
    lenomega <- length(omegahat)
    leneta <- getleneta(cov.model)
    lenY <- length(Yhat)
    npars <- lenbeta + lenomega + leneta + lenY
    
    sigma <- matrix(0,npars,npars)
    
    # eta
    logpost <- function(eta,surv,X,beta,omega,Y,priors,cov.model,u,control){

        sigma <- matrix(EvalCov(cov.model=cov.model,u=u,parameters=eta),n,n)
        cholsigma <- t(chol(sigma))
        cholsigmainv <- solve(cholsigma)
        MU <- -cov.model$itrans[[control$sigmaidx]](eta[control$sigmaidx])^2/2
        gamma <- cholsigmainv%*%(Y-MU)  
        
        logpost <- logPosterior(surv=surv,X=X,beta=beta,omega=omega,eta=eta,gamma=gamma,priors=priors,cov.model=cov.model,u=u,control=control)        
                          
        return(logpost)
    }

    npts <- 20
    if(leneta>=3){
        npts <- 10
    }
    rgs <- getparranges(priors=priors,leneta=leneta)   
    qa <- QuadApprox(logpost,npts=npts,argRanges=rgs,surv=surv,X=X,beta=betahat,omega=omegahat,Y=Yhat,priors=priors,cov.model=cov.model,u=u,control=control)
    
    matr <- qa$curvature
    etahat <- qa$max
    
    # entry for eta in propossal covariance
    sigma[(lenbeta+lenomega+1):(lenbeta+lenomega+leneta),(lenbeta+lenomega+1):(lenbeta+lenomega+leneta)] <- matr    
        
    #estimate of gamma  
    ssigma <- matrix(EvalCov(cov.model=cov.model,u=u,parameters=etahat),n,n)
    cholssigma <- t(chol(ssigma))
    MU <- -cov.model$itrans[[control$sigmaidx]](etahat[control$sigmaidx])^2/2
    gammahat <- solve(cholssigma)%*%(Yhat-MU)
    
    hessian <- logPosterior(surv=surv,X=X,beta=betahat,omega=omegahat,eta=etahat,gamma=gammahat,priors=priors,cov.model=cov.model,u=u,control=control,hessian=TRUE)
    
    # beta and omega
    sigma[1:lenbeta,1:lenbeta] <- hessian$hess_beta
    sigma[(lenbeta+1):(lenbeta+lenomega),(lenbeta+1):(lenbeta+lenomega)] <- hessian$hess_omega
    sigma[(lenbeta+1):(lenbeta+lenomega),(1:lenbeta)] <- hessian$hess_omega_beta
    sigma[(1:lenbeta),(lenbeta+1):(lenbeta+lenomega)] <- t(hessian$hess_omega_beta)       
    # gamma
    diag(sigma)[(lenbeta+lenomega+leneta+1):npars] <- hessian$hess_gamma   
    
    return(list(etahat=etahat,sigma=solve(-sigma))) 
}




##' proposalVariance_gridded function
##'
##' A function to 
##'
##' @param X X 
##' @param surv X 
##' @param betahat X 
##' @param omegahat X 
##' @param Yhat X 
##' @param priors X 
##' @param cov.model X 
##' @param u X 
##' @param control X 
##' @return ...
##' @export

proposalVariance_gridded <- function(X,surv,betahat,omegahat,Yhat,priors,cov.model,u,control){

    Ygrid <- gridY(Y=Yhat,control=control)    
     
    n <- nrow(X)
    lenbeta <- length(betahat)
    lenomega <- length(omegahat)
    leneta <- getleneta(cov.model)
    lenY <- length(Ygrid)
    npars <- lenbeta + lenomega + leneta + lenY
    
    # eta
    logpost <- function(eta,surv,X,beta,omega,Ygrid,priors,cov.model,u,control){
        
        covbase <- matrix(EvalCov(cov.model=cov.model,u=u,parameters=eta),control$Mext,control$Next)
        
        rootQeigs <- sqrt(1/Re(fft(covbase)))   
       
        pars <- sapply(1:length(eta),function(i){cov.model$itrans[[i]](eta[i])})
        ymean <- -pars[which(cov.model$parnames=="sigma")]^2/2
        gamma <- GammafromY(Ygrid,rootQeigs=rootQeigs,mu=ymean)  
        
        logpost <- logPosterior_gridded(surv=surv,X=X,beta=beta,omega=omega,eta=eta,gamma=gamma,priors=priors,cov.model=cov.model,u=u,control=control)        
                          
        return(logpost)
    }

    npts <- 20
    if(leneta>=3){
        npts <- 10
    }
    rgs <- getparranges(priors=priors,leneta=leneta)   
    qa <- QuadApprox(logpost,npts=npts,argRanges=rgs,surv=surv,X=X,beta=betahat,omega=omegahat,Ygrid=Ygrid,priors=priors,cov.model=cov.model,u=u,control=control)
    
    matr <- qa$curvature
    etahat <- qa$max
    
    sigma <- matrix(0,lenbeta + lenomega + leneta,lenbeta + lenomega + leneta)
    
    # entry for eta in propossal covariance
    sigma[(lenbeta+lenomega+1):(lenbeta+lenomega+leneta),(lenbeta+lenomega+1):(lenbeta+lenomega+leneta)] <- matr    
        
    #estimate of gamma  
    covbase <- matrix(EvalCov(cov.model=cov.model,u=u,parameters=etahat),control$Mext,control$Next)        
    rootQeigs <- sqrt(1/Re(fft(covbase)))   
    pars <- sapply(1:length(etahat),function(i){cov.model$itrans[[i]](etahat[i])})
    ymean <- -pars[which(cov.model$parnames=="sigma")]^2/2
    gammahat <- GammafromY(Ygrid,rootQeigs=rootQeigs,mu=ymean)
    
    hessian <- logPosterior_gridded(surv=surv,X=X,beta=betahat,omega=omegahat,eta=etahat,gamma=gammahat,priors=priors,cov.model=cov.model,u=u,control=control,hessian=TRUE)
    
    # beta and omega
    sigma[1:lenbeta,1:lenbeta] <- hessian$hess_beta
    sigma[(lenbeta+1):(lenbeta+lenomega),(lenbeta+1):(lenbeta+lenomega)] <- hessian$hess_omega
    sigma[(lenbeta+1):(lenbeta+lenomega),(1:lenbeta)] <- hessian$hess_omega_beta
    sigma[(1:lenbeta),(lenbeta+1):(lenbeta+lenomega)] <- t(hessian$hess_omega_beta)       
    # gamma
    hess_gam <- hessian$hess_gamma 
    
    sigma <- (-1) * sigma # variance is inverse of observed information    
    
    matidx <- (lenbeta+lenomega+leneta+1):npars
    matidx <- matrix(matidx,nrow=length(matidx),ncol=2) 

    sigmaret <- Matrix(0,npars,npars)
    sigmaret[1:(lenbeta+lenomega+leneta),1:(lenbeta+lenomega+leneta)] <- solve(sigma)
    sigmaret[matidx] <- -1/hess_gam   
    
    return(list(etahat=etahat,sigma=sigmaret)) 
}



##' estimateY function
##'
##' A function to 
##'
##' @param X X 
##' @param betahat X 
##' @param omegahat X 
##' @param surv X
##' @param control X 
##' @return ...
##' @export

estimateY <- function(X,betahat,omegahat,surv,control){

    censoringtype <- attr(surv,"type")
   
    omega <- control$omegaitrans(omegahat) # this is omega on the correct scale
    
    haz <- setupHazard(dist=control$dist,pars=omega,grad=FALSE,hess=FALSE)

    n <- nrow(X)
    
    if(censoringtype=="left" | censoringtype=="right"){
        notcensored <- surv[,"status"]==1
    }
    else{
        rightcensored <- surv[,"status"] == 0
        notcensored <- surv[,"status"] == 1
        lefttruncated <- surv[,"status"] == 2
        intervalcensored <- surv[,"status"] == 3
    }   

    # setup function J=exp(X%*%beta + Y)*H_0(t)
    if(censoringtype=="left" | censoringtype=="right"){
        tsubs <- surv[,"time"]
    }
    else{ # else interval censored
        tsubs <- surv[,"time1"]        
    } 
    
    for(i in 1:n){
    
        if(notcensored[i]){
            next
        }

        if(censoringtype=="left" | censoringtype=="right"){            
            if(censoringtype=="right"){
                tpot <- tsubs[notcensored][tsubs[notcensored]>tsubs[i]] # potential t                    
            }
            else{ # censoringtype=="left"
                tpot <- tsubs[notcensored][tsubs[notcensored]<tsubs[i]] # potential t
            }
        }
        else{
            if(rightcensored[i]){
                tpot <- tsubs[notcensored][tsubs[notcensored]>tsubs[i]] # potential t
            }
            if(lefttruncated[i]){
                tpot <- tsubs[notcensored][tsubs[notcensored]<tsubs[i]] # potential t
            }
            if(intervalcensored[i]){
                tpot <- surv[,"time1"][i] + 0.5*(surv[,"time2"][i]-surv[,"time1"][i]) # mid point of interval
            }
        }
            
        if(length(tpot)==0){
            next # leave tsubs[i] alone 
        }
        else{
            tsubs[i] <- sample(tpot,1) # ignoring covariates, sample from empirical distribution of times exceeding (right censored), or less than (left censored) the observed time 
        }

    }
    
    
    
    Y <- -X%*%betahat - log(haz$H(tsubs)) # greedy estimate of Y (maximise individual contributions to log-likelihood) ... note log(delta) is now omitted  

    return(Y)    
}