##' getcov function
##'
##' A function to return the covariance from a model based on the randomFields covariance functions. Not intended for general use.
##'
##' @param u distance
##' @param sigma variance parameter
##' @param phi scale parameter
##' @param model correlation type, see ?CovarianceFct
##' @param pars vector of additional parameters for certain classes of covariance function (eg Matern), these must be supplied in the order given in ?CovarianceFct and are not estimated
##' @return this is just a wrapper for CovarianceFct
##' @export

getcov <- function(u,sigma,phi,model,pars){
    return(suppressWarnings(CovarianceFct(x=u,param=c(mean=0,variance=sigma^2,nugget=0,scale=phi,pars),model=model)))
}



##' covmodel function
##'
##' A function to define the spatial covariance model, see also ?CovarianceFct. Note that the parameters defined by the 'pars' argument are fixed,
##' i.e. not estimated by the MCMC algorithm. To have spatsurv estimate these parameters, the user must construct a new covariance function to do so, 
##' see the spatsurv vignette. 
##'
##' @param model correlation type, a string see ?CovarianceFct 
##' @param pars vector of additional parameters for certain classes of covariance function (eg Matern), these must be supplied in the order given in ?CovarianceFct and are not estimated 
##' @return an object of class covmodel
##' @seealso \link{CovarianceFct}
##' @export

covmodel <- function(model,pars){
    retlist <- list()
    retlist$model <- model
    retlist$pars <- pars
    retlist$npar <- 2
    retlist$parnames <- c("sigma","phi")
    retlist$itrans <- exp # inverse transform back to correct scale
    retlist$trans <- log # transform assumed   
    class(retlist) <- c("covmodel","fromRandomFieldsCovarianceFct")
    return(retlist)
}


##' ExponentialCovFct function
##'
##' A function to declare and also evaluate an exponential covariance function.
##'
##' @return the exponential covariance function
##' @seealso \link{SpikedExponentialCovFct}, \link{covmodel}, \link{CovarianceFct}
##' @export

ExponentialCovFct <- function(){
    ans <- list()
    ans$npar <- 2
    ans$parnames <- c("sigma","phi")
    ans$itrans <- exp # inverse transform back to correct scale
    ans$trans <- log # transform assumed  
    ans$eval <- function(u,pars){
        ans<- pars[1]^2 * exp(-u/pars[2])        
        return(ans)
    }
    class(ans) <- c("covmodel","fromUserFunction")
    return(ans)
}

##' SpikedExponentialCovFct function
##'
##' A function to declare and also evaluate a spiked exponential covariance function. This is an exponential covariance function with a nugget.
##'
##' @return the spiked exponential covariance function
##' @seealso \link{ExponentialCovFct}, \link{covmodel}, \link{CovarianceFct}
##' @export

SpikedExponentialCovFct <- function(){
    ans <- list()
    ans$npar <- 3
    ans$parnames <- c("sigma","phi","nugget")
    ans$itrans <- exp # inverse transform back to correct scale
    ans$trans <- log # transform assumed  
    ans$eval <- function(u,pars){
        ans <- pars[1]^2 * exp(-u/pars[2])
        ans[u==0] <- ans[u==0] + pars[3]^2 
        return(ans)
    }
    class(ans) <- c("covmodel","fromUserFunction")
    return(ans)
}




##' EvalCov function
##'
##' This function is used to evaluate the covariance function within the MCMC run. Not intended for general use.
##'
##' @param cov.model an object of class covmodel
##' @param u vector of distances
##' @param parameters vector of parameters
##' @return method EvalCov
##' @export

EvalCov <- function(cov.model,u,parameters){
    if(inherits(cov.model,"fromRandomFieldsCovarianceFct")){
        ev <- getcov(u=u,sigma=parameters[1],phi=parameters[2],model=cov.model$model,pars=cov.model$pars)
    }
    else if(inherits(cov.model,"fromUserFunction")){
        ev <- cov.model$eval(u=u,pars=parameters)
    }
    else{
        stop("Unknown covariance type")
    }
    return(ev)
}


