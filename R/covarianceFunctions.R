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
    if(model=="exponential"){
        return(sigma^2*exp(-u/phi))
    }
    else if(model=="Matern32"){
        return(sigma^2*(1+sqrt(3)*u/phi)*exp(-sqrt(3)*u/phi))
    }
    else if(model=="Matern52"){
        return(sigma^2*(1+sqrt(5)*u/phi+5*u^2/(3*phi^2))*exp(-sqrt(5)*u/phi))
    }
    else{
        stop("Functionality is temporarily unavailable due to deprecation of RandomFields package. Models: 'exponential', 'Matern32' and 'Matern52' are currently available")
    }
    #return(suppressWarnings(CovarianceFct(x=u,param=c(mean=0,variance=sigma^2,nugget=0,scale=phi,pars),model=model)))
}



##' covmodel function
##'
##' A function to define the spatial covariance model, see also ?CovarianceFct. Note that the parameters defined by the 'pars' argument are fixed,
##' i.e. not estimated by the MCMC algorithm. To have spatsurv estimate these parameters, the user must construct a new covariance function to do so, stop("")
##' see the spatsurv vignette.
##'
##' @param model correlation type, a string see ?CovarianceFct
##' @param pars vector of additional parameters for certain classes of covariance function (eg Matern), these must be supplied in the order given in ?CovarianceFct and are not estimated
##' @return an object of class covmodel
## @seealso CovarianceFct
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
##' @seealso \link{SpikedExponentialCovFct}, \link{covmodel}
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

##' Independent function
##'
##' A function to declare and also evaluate an exponential covariance function.
##'
##' @return the exponential covariance function
##' @seealso \link{SpikedExponentialCovFct}, \link{covmodel}
##' @export

Independent <- function(){
    ans <- list()
    ans$npar <- 1
    ans$parnames <- "sigma"
    ans$itrans <- exp # inverse transform back to correct scale
    ans$trans <- log # transform assumed
    ans$eval <- function(u,pars){
        ans<- rep(pars[1]^2,length(u))
        ans[u>0] <- 0
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
##' @seealso \link{ExponentialCovFct}, \link{covmodel}
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


##' SPDE function
##'
##' A function to declare and evaluate an SPDE covariance function.
##'
##' @param ord the order of the model to be used, currently an integer between 1 an 3. See Lindgren 2011 paper.
##' @return an covariance function based on the SPDE model
##' @seealso \link{ExponentialCovFct}, \link{covmodel}
##' @export

SPDE <- function(ord){
    ans <- list()
    ans$npar <- 2
    ans$parnames <- c("psi","a")
    ans$itrans <- function(x){return(c(exp(x[1]),4+exp(x[2])))} # inverse transform back to correct scale
    ans$trans <- function(x){return(c(log(x[1]),log(x[2]-4)))} # transform assumed
    ans$eval <- function(precMatStruct,pars){
        ans <- list()
        ans$precmat <- precMatStruct(SPDEprec(pars[2],ord))
        ans$vmult <- pars[1]
        return(ans)
    }
    ans$order <- ord
    class(ans) <- c("covmodel","SPDEmodel")
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
