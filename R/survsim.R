##' simsurv function
##'
##' A function to simulate spatial parametric proportional hazards model with baseline hazard derived from the exponential or weibull model. The function works
##' by simulating candidate survival times using MCMC in parallel for each individual based on each individual's covariates and the common
##' parameter effects, beta.  
##'
##' @param X a matrix of covariate information 
##' @param beta the parameter effects 
##' @param theta parameter for the baseline hazard model (the rate for exponential data and the shape and scale for Weibull data)
##' @param dist the distribution choice: exp or weibull at present
##' @param coords matrix with 2 columns giving the coordinates at which to simulate data
##' @param sigmaphi a vector of length 2: the parameters for the covariance function, sigma and phi in that order
##' @param cov.model an object of class covmodel, see ?covmodel
##' @param mcmc.control mcmc control paramters, see ?mcmcpars
##' @param savechains save all chains? runs faster if set to FALSE, but then you'll be unable to conduct convergence/mixing diagnostics
##' @return simulated survival times from the exponential model (the last simulated value from the MCMC chains)
##' @seealso \link{covmodel}, \link{survspat} 
##' @export

simsurv <- function(X=cbind(age=runif(100,5,50),sex=rbinom(100,1,0.5),cancer=rbinom(100,1,0.2)),
                            beta=matrix(c(0.0296,0.0261,0.035),3,1),
                            theta=1,
                            dist="exp",
                            coords=matrix(runif(2*nrow(X)),nrow(X),2),
                            sigmaphi=c(1,0.1),
                            cov.model=covmodel(model="exponential",pars=NULL),
                            mcmc.control=mcmcpars(nits=100000,burn=10000,thin=90),      
                            savechains=TRUE){

    mcmcloop <- mcmcLoop(N=mcmc.control$nits,burnin=mcmc.control$burn,thin=mcmc.control$thin,progressor=mcmcProgressTextBar)
    
    distmat <- as.matrix(dist(coords))

    n <- nrow(X)
    u <- as.vector(distmat)
    
    sigma <- matrix(getcov(u=u,sigma=sigmaphi[1],phi=sigmaphi[2],model=cov.model$model,pars=cov.model$pars),n,n)
    sigmachol <- t(chol(sigma))
    Y <- -sigmaphi[1]^2/2+sigmachol%*%rnorm(n)
    expY <- exp(Y)  
 
    XbetaplusY <- X%*%beta + Y
    expXbetaplusY <- exp(XbetaplusY)
    #rexpXbeta <- rate*expXbeta
    
    nmatrows <- ceiling((mcmc.control$nits-mcmc.control$burn-mcmcloop$waste)/mcmc.control$thin)
    
    if(savechains){
        T <- matrix(NA,nmatrows,n) #matrix to store simulated times
    }
    else{
        T <- NULL
    }
    tarrec <- matrix(NA,nmatrows,n)
    acrec <- rep(NA,nmatrows)
    count <- 1 # counter to index retained iteration numbers
    
    if(dist=="exp"){
        t <- rexp(n,theta)
    }
    else if(dist=="weibull"){
        t <- rweibull(n,shape=theta[1],scale=theta[2])
    }
    
    oldltar <- do.call(paste(dist,"_ltar",sep=""),list(t=t,XbetaplusY=XbetaplusY,expXbetaplusY=expXbetaplusY,theta=theta))
    
    while(nextStep(mcmcloop)){        
    
        newt <- rexp(n,rate=1/t) 
        newltar <- do.call(paste(dist,"_ltar",sep=""),list(t=newt,XbetaplusY=XbetaplusY,expXbetaplusY=expXbetaplusY,theta=theta))
        
        frac <- exp(newltar-oldltar+dexp(t,1/newt)-dexp(newt,1/t))        
        
        ac <- pmin(1,frac)
        
        keepnew <- ac>runif(n)
        t[keepnew] <- newt[keepnew]
        oldltar[keepnew] <- newltar[keepnew]
    
        if(is.retain(mcmcloop)){
            if(savechains){        
                T[count,] <- t
            }
            tarrec[count,] <- oldltar
            acrec[count] <- mean(ac)
            count <- count + 1
        }        
    }   
    
    if(is.na(tarrec[nmatrows])){
        if(savechains){
            T <- T[-nmatrows,]
        }
        tarrec <- tarrec[-nmatrows]
        acrec <- acrec[-nmatrows]
    }
    
    if(savechains){
        tchoice <- T[nrow(T),]
    }
    else{
        tchoice <- t
    }
    
    cat("Mean acceptance:",mean(acrec),"\n")
    
    return(list(X=X,T=T,survtimes=tchoice,theta=theta,beta=beta,tarrec=tarrec,n=n,Y=Y,dist=dist,coords=coords,sigmaphi=sigmaphi,distmat=distmat,u=u))
}