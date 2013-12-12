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



##
## NOTE THIS FUNCTION HAS BEEN REPLACED WITH QuadApprox
##
## quadapprox function
##
## NOTE THIS FUNCTION HAS NOW BEEN SUPERCEDED BY THE FUNCTION QuadApprox
##
##
## A function to compute the second derivative of a function using a quadratic approximation to the function on a 
## grid of points defined by xseq and yseq. Also returns the local maximum. 
##
## @param fun a function 
## @param xseq sequence of x-values defining grid on which to compute the approximation 
## @param yseq sequence of y-values defining grid on which to compute the approximation 
## @param ... other arguments to be passed to fun
## @return a 2 by 2 matrix containing the curvature at the maximum and the (x,y) value at which the maximum occurs 
## @export


#quadapprox <- function(fun,xseq,yseq,...){
#    nx <- length(xseq)
#    ny <- length(yseq)
#    funvals <- matrix(NA,nx,ny)
#    #pb <- txtProgressBar(min=0,max=nx*ny,style=3)
#    count <- 0
#    gr <- expand.grid(xseq,yseq)
#    funvals <- apply(gr,1,function(xy){fun(c(xy[1],xy[2]),...)})
#    image.plot(xseq,yseq,matrix(funvals,nx,ny),main="Approx Posterior")
#    funvals <- as.vector(funvals)
#    x <- gr[,1]
#    x2 <- gr[,1]^2
#    y <- gr[,2]    
#    y2 <- gr[,2]^2
#    xy <- gr[,1]*gr[,2]
#    mod <- lm(funvals~x2+x+y2+y+xy)
#    image.plot(xseq,yseq,matrix(fitted(mod),nx,ny),main="Quadratic Approximation")
#    co <- coefficients(mod)
#    d2dxx <- 2*co[2]
#    d2dyy <- 2*co[4]
#    d2dydx <- co[6]   
#    
#    sigmainv <- matrix(c(d2dxx,d2dydx,d2dydx,d2dyy),2,2) 
#    
#    #print(-solve(sigmainv))
#    
#    eta1est <- (-2*co[3]*co[4]+co[5]*co[6])/(4*co[2]*co[4]-co[6]^2)
#    eta2est <- (-co[5]-co[6]*eta1est)/(2*co[4])
#    
#    #print(paste("Estimated exp(eta):",exp(eta1est),",",exp(eta2est)))
#    
#    return(list(max=c(eta1est,eta2est),curvature=sigmainv,mod=mod)) 
#}



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