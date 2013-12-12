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




##' quadapprox function
##'
##' A function to compute the second derivative of a function using a quadratic approximation to the function on a 
##' grid of points defined by xseq and yseq. Also returns the local maximum. 
##'
##' @param fun a function 
##' @param xseq sequence of x-values defining grid on which to compute the approximation 
##' @param yseq sequence of y-values defining grid on which to compute the approximation 
##' @param ... other arguments to be passed to fun
##' @return a 2 by 2 matrix containing the curvature at the maximum and the (x,y) value at which the maximum occurs 
##' @export


quadapprox <- function(fun,xseq,yseq,...){
    nx <- length(xseq)
    ny <- length(yseq)
    funvals <- matrix(NA,nx,ny)
    #pb <- txtProgressBar(min=0,max=nx*ny,style=3)
    count <- 0
    gr <- expand.grid(xseq,yseq)
    funvals <- apply(gr,1,function(xy){fun(c(xy[1],xy[2]),...)})
    image.plot(xseq,yseq,matrix(funvals,nx,ny),main="Approx Posterior")
    funvals <- as.vector(funvals)
    x <- gr[,1]
    x2 <- gr[,1]^2
    y <- gr[,2]    
    y2 <- gr[,2]^2
    xy <- gr[,1]*gr[,2]
    mod <- lm(funvals~x2+x+y2+y+xy)
    image.plot(xseq,yseq,matrix(fitted(mod),nx,ny),main="Quadratic Approximation")
    co <- coefficients(mod)
    d2dxx <- 2*co[2]
    d2dyy <- 2*co[4]
    d2dydx <- co[6]   
    
    sigmainv <- matrix(c(d2dxx,d2dydx,d2dydx,d2dyy),2,2) 
    
    print(-solve(sigmainv))
    
    eta1est <- (-2*co[3]*co[4]+co[5]*co[6])/(4*co[2]*co[4]-co[6]^2)
    eta2est <- (-co[5]-co[6]*eta1est)/(2*co[4])
    
    print(paste("Estimated exp(eta):",exp(eta1est),",",exp(eta2est)))
    
    return(list(max=c(eta1est,eta2est),curvature=sigmainv,mod=mod)) 
}