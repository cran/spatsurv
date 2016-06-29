##' imputationModel function
##'
##' A function to 
##'
##' @param formula X 
##' @param offset X 
##' @param covariateData X 
##' @param priors X 
##' @return ...
##' @export

imputationModel <- function(formula,offset,covariateData,priors){
    retlist <- list()
    retlist$formula <- formula
    retlist$offset <- offset
    retlist$covariateData <- covariateData
    retlist$priors <- priors
    return(retlist)
}



