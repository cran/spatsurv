##' getcov function
##' @param u distance
##' @param sigma variance parameter
##' @param phi scale parameter
##' @param model correlation type, see ?CovarianceFct
##' @param pars vector of additional parameters for certain classes of covariance function (eg Matern), these must be supplied in the order given in ?CovarianceFct and are not estimated
##' @return this is just a wrapper for CovarianceFct

getcov <- function(u,sigma,phi,model,pars){
    return(suppressWarnings(CovarianceFct(x=u,param=c(mean=0,variance=sigma^2,nugget=0,scale=phi,pars),model=model)))
}



##' covmodel function
##'
##' A function to define the spatial covariance model, see also ?getcov and ?CovarianceFct
##'
##' @param model correlation type, see ?CovarianceFct 
##' @param pars vector of additional parameters for certain classes of covariance function (eg Matern), these must be supplied in the order given in ?CovarianceFct and are not estimated 
##' @return an object of class covmodel
##' @export

covmodel <- function(model,pars){
    retlist <- list()
    retlist$model <- model
    retlist$pars <- pars
    retlist$npar <- 2
    retlist$parnames <- c("sigma","phi")
    retlist$itrans <- c(exp,exp) # inverse transform back to correct scale
    retlist$trans <- c(log,log) # transform assumed   
    class(retlist) <- c("covmodel","fromRandomFieldsCovarianceFct")
    return(retlist)
}


##' ExponentialCovFct function
##'
##' A function to declare and also evaluate an exponential covariance function.
##'
## @param u distance
## @param pars parameters of the latent field, a vector. The first parameter is the marginal standard deviation i.e. sigma.
##' @return the exponential covariance function
##' @export

ExponentialCovFct <- function(){
    ans <- list()
    ans$npar <- 2
    ans$parnames <- c("sigma","phi")
    ans$itrans <- c(exp,exp) # inverse transform back to correct scale
    ans$trans <- c(log,log) # transform assumed  
    ans$eval <- function(u,pars){
        pars <- exp(pars)
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
## @param u distance
## @param pars parameters of the latent field, a vector. The first parameter is the marginal standard deviation i.e. sigma.
##' @return the spiked exponential covariance function
##' @export

SpikedExponentialCovFct <- function(){
    ans <- list()
    ans$npar <- 3
    ans$parnames <- c("sigma","phi","nugget")
    ans$itrans <- c(exp,exp,exp) # inverse transform back to correct scale
    ans$trans <- c(log,log,log) # transform assumed  
    ans$eval <- function(u,pars){
        pars <- exp(pars)
        ans <- pars[1]^2 * exp(-u/pars[2])
        ans[u==0] <- ans[u==0] + pars[3]^2 
        return(ans)
    }
    class(ans) <- c("covmodel","fromUserFunction")
    return(ans)
}




##' EvalCov function
##'
##' Evaluate the covariance function
##'
##' @param cov.model an object of class covmodel
##' @param u vector of distances
##' @param parameters vector of parameters
##' @return method EvalCov
##' @export

EvalCov <- function(cov.model,u,parameters){
    if(inherits(cov.model,"fromRandomFieldsCovarianceFct")){
        ev <- getcov(u=u,sigma=exp(parameters[1]),phi=exp(parameters[2]),model=cov.model$model,pars=cov.model$pars)
    }
    else if(inherits(cov.model,"fromUserFunction")){
        ev <- cov.model$eval(u=u,pars=parameters)
    }
    else{
        stop("Unknkown covariance type")
    }
    return(ev)
}


