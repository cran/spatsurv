##' gensens function
##'
##' A function to generate observed times given a vector of true survival times and a vector of censoring times. Used in the simulation of
##' survival data
##'
##' @param survtimes a vector of survival times 
##' @param censtimes a vector of censoring times 
##' @return a named list containing 'obstimes', the observed time of the event; and 'censored', the censoring indicator which is equal to 1 if the
##' event is observed and 0 otherwise.
##' @export


gensens <- function(survtimes,censtimes){
    
    n <- length(survtimes)
    
    if(length(survtimes)!=length(censtimes)){
        stop("survtimes and censtimes should have the same length")
    } 

    obstimes <- survtimes
    cens <- rep(1,n)
   
    for(i in 1:n){
        if(censtimes[i]<survtimes[i]){
            obstimes[i] <- censtimes[i]
            cens[i] <- 0
        }
    }
    
    return(Surv(time=obstimes,event=cens))
}



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
    class(retlist) <- "covmodel"
    return(retlist)
}



#empvari <- function(coords,z,plot=TRUE,...){
#    if(length(z)!=nrow(coords)){
#        stop("length(z) must equal nrow(coords)")
#    }
#    distmat <- as.matrix(dist(coords))
#    diffz2 <- outer(as.vector(z),as.vector(z),FUN="-")^2 # squared differences between z values
#    lt <- lower.tri(distmat)
#
#    x <- as.vector(distmat[lt])
#    y <- as.vector(diffz2[lt])
#    
#    if(plot){
#        plot(x,y,pch=".")  
#        lines(lowess(x,y,...),col="red")
#    } 
#    
#    browser()
#    return(list(x=x,y=y))
#}