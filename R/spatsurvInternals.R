##' getleneta function
##'
##' A function to compute the length of eta
##'
##' @param cov.model a covariance model 
##' @return the length of eta

getleneta <- function(cov.model){
    if(inherits(cov.model,"fromRandomFieldsCovarianceFct")){
        leneta <- 2
    }
    else if(inherits(cov.model,"fromUserFunction")){
        leneta <- cov.model$npar
    }
    else{
        stop("Unknkown covariance type")
    }
    return(leneta)
} 


##' getparranges function
##'
##' A function to 
##'
##' @param priors X 
##' @param leneta X 
##' @param mult X
##' @return ...
##' @export

getparranges <- function(priors,leneta,mult=1.96){
    rgs <- list()
    if(length(priors$etaprior$mean)==1){
        for(i in 1:leneta){
            rgs[[i]] <- c(priors$etaprior$mean-mult*priors$etaprior$sd,priors$etaprior$mean+mult*priors$etaprior$sd)
        }
    }
    else{
        for(i in 1:leneta){
            rgs[[i]] <- c(priors$etaprior$mean[i]-mult*priors$etaprior$sd[i],priors$etaprior$mean[i]+mult*priors$etaprior$sd[i])
        }
    } 
    return(rgs)
}    