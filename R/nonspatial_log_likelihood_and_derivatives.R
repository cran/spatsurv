##' maxlikparamPHsurv function
##'
##' A function to get initial estimates of model parameters using maximum likelihood. Not intended for general purose use.
##'
##' @param surv an object of class Surv
##' @param X the design matrix, containing covariate information 
##' @param control a list containg various control parameters for the MCMC and post-processing routines  
##' @return initial estimates of the parameters
##' @export

maxlikparamPHsurv <- function(surv,X,control){
    
    
    likfun <- function(pars){
        beta <- pars[1:ncol(X)]
        omega <- pars[(ncol(X)+1):length(pars)]
        return(NonSpatialLogLikelihood_or_gradient(surv=surv,X=X,beta=beta,omega=omega,control=control,loglikelihood=TRUE,gradient=FALSE))
    }
    
    gradfun <- function(pars){
        beta <- pars[1:ncol(X)]
        omega <- pars[(ncol(X)+1):length(pars)]
        return(NonSpatialLogLikelihood_or_gradient(surv=surv,X=X,beta=beta,omega=omega,control=control,loglikelihood=FALSE,gradient=TRUE))
    }
    
    betainit <- rep(0,ncol(X))
    
    if(is.null(distinfo(control$dist)()$MLinits)){
        omegainit <- rep(1e-10,distinfo(control$dist)()$npars)
    }
    else{
        omegainit <- distinfo(control$dist)()$MLinits
    }
    
    #browser()
    opt <- try(optim(par=c(betainit,omegainit),fn=likfun,gr=gradfun,method="BFGS",control=control$optimcontrol,hessian=control$hessian))
    if(inherits(opt,"try-error")){
        stop("Possible problem with initial values in obtaining maximum likelihood estimates of parameters, try setting MLinits in function distinfo for chosen baseline hazard distribution")
    }
    
    return(opt)
}


##' NonSpatialLogLikelihood_or_gradient function
##'
##' A function to evaluate the log-likelihood of a non-spatial parametric proportional hazards model. Not intended for general use.
##'
##' @param surv an object of class Surv 
##' @param X the design matrix, containing covariate information 
##' @param beta parameter beta 
##' @param omega parameter omega 
##' @param control a list containg various control parameters for the MCMC and post-processing routines   
##' @param loglikelihood logical whether to evaluate the log-likelihood
##' @param gradient logical whether to evaluate the gradient
##' @return ...
##' @export

NonSpatialLogLikelihood_or_gradient <- function(surv,X,beta,omega,control,loglikelihood,gradient){
    
    censoringtype <- attr(surv,"type")
    
    omegaorig <- omega # recall we are working with omega on the transformed scale
    omega <- control$omegaitrans(omega) # this is omega on the correct scale
    
    haz <- setupHazard(dist=control$dist,pars=omega,grad=TRUE)
    
    n <- nrow(X)
    
    Xbeta <- X%*%beta
    expXbeta <- exp(Xbeta)
    
    if(censoringtype=="left" | censoringtype=="right"){
        censored <- surv[,"status"]==0
        notcensored <- !censored

        Ctest <- any(censored)
        Utest <- any(notcensored)        
        
    }
    else{
        rightcensored <- surv[,"status"] == 0
        notcensored <- surv[,"status"] == 1
        leftcensored <- surv[,"status"] == 2
        intervalcensored <- surv[,"status"] == 3

        Rtest <- any(rightcensored)        
        Utest <- any(notcensored) 
        Ltest <- any(leftcensored)
        Itest <- any(intervalcensored)
    }

    

    # setup function J=exp(X%*%beta)*H_0(t)
    if(censoringtype=="left" | censoringtype=="right"){
        H <- haz$H(surv[,"time"])
        h <- haz$h(surv[,"time"])
        J <- expXbeta*H
        
        J <- as.vector(J)
    }
    else{ # else interval censored
        H1 <- haz$H(surv[,"time1"])
        H2 <- haz$H(surv[,"time2"])
        J1 <- expXbeta*H1
        J2 <- expXbeta*H2
        
        J1 <- as.vector(J1)
        J2 <- as.vector(J2)
    }    
    
    # setup function S=exp(-J(t))
    if(censoringtype=="right"|censoringtype=="left"){ 
        S <- exp(-J)
    }    
    else{
        S1 <- exp(-J1)
        S2 <- exp(-J2)
    }
    
    if(censoringtype=="right" | censoringtype=="left"){
        h <- haz$h(surv[,"time"])
    }
    else{ #censoringtype=="interval"
        h <- haz$h(surv[,"time1"])
    }
    
    
    if(loglikelihood){ 
        if(censoringtype=="right"){
            loglik <-  (if(Utest){sum(Xbeta[notcensored] + log(h)[notcensored] - J[notcensored])}else{0}) + 
                       (if(Ctest){sum(-J[censored])}else{0})
        }
        else if(censoringtype=="left"){
            loglik <-  (if(Utest){sum(Xbeta[notcensored] + log(h)[notcensored] - J[notcensored])}else{0}) + 
                       (if(Ctest){sum(log(1-S[censored]))}else{0})
        }
        else{ #censoringtype=="interval"

            loglik <-  (if(Utest){sum(Xbeta[notcensored] + log(h)[notcensored] - J1[notcensored])}else{0}) + 
                       (if(Rtest){sum(-J1[rightcensored])}else{0}) + 
                       (if(Ltest){sum(log(1-S1[leftcensored]))}else{0}) + 
                       (if(Itest){sum(log(S1[intervalcensored]-S2[intervalcensored]))}else{0})
        }
        
        return(-loglik)
    }
    
    if(gradient){
    
        if(censoringtype=="left" | censoringtype=="right"){
            dJ_dbeta <- J*X
            dJ_domega <- matrix(as.vector(expXbeta)*haz$gradH(surv[,"time"]),ncol=length(omega))
        }
        else{ # interval censoring
            dJ_dbeta1 <- J1*X
            dJ_domega1 <- matrix(as.vector(expXbeta)*haz$gradH(surv[,"time1"]),ncol=length(omega))
            
            dJ_dbeta2 <- J1*X
            dJ_domega2 <- matrix(as.vector(expXbeta)*haz$gradH(surv[,"time2"]),ncol=length(omega))
        }        
            
        if(censoringtype=="right"){
            dP_dbeta <- (if(Utest){colSums(X[notcensored,,drop=FALSE] - dJ_dbeta[notcensored,,drop=FALSE])}else{0}) + 
                        (if(Ctest){colSums(-dJ_dbeta[censored,,drop=FALSE])}else{0})
            dh_domega <- matrix(haz$gradh(surv[,"time"]),ncol=length(omega))
            dP_domega <-    (if(Utest){colSums(dh_domega[notcensored,,drop=FALSE] / h[notcensored] - dJ_domega[notcensored,])}else{0}) + 
                            (if(Ctest){colSums(-dJ_domega[censored,,drop=FALSE])}else{0})
            dP_domega <- control$omegajacobian(omegaorig)*dP_domega # this puts the derivative back on the correct scale dL/dpsi = dL/dtheta * dtheta/dpsi, e.g. psi=log(theta)

            grad <- c(dP_dbeta,dP_domega)
            
        }
        else if(censoringtype=="left"){
            dP_dbeta <- (if(Utest){colSums(X[notcensored,,drop=FALSE] - dJ_dbeta[notcensored,,drop=FALSE])}else{0}) + 
                        (if(Ctest){colSums((S[censored]/(1-S[censored]))*dJ_dbeta[censored,,drop=FALSE])}else{0})
            dh_domega <- matrix(haz$gradh(surv[,"time"]),ncol=length(omega))
            dP_domega <-    (if(Utest){colSums(dh_domega[notcensored,,drop=FALSE] / h[notcensored] - dJ_domega[notcensored,,drop=FALSE])}else{0}) + 
                            (if(Ctest){colSums((S[censored]/(1-S[censored]))*dJ_domega[censored,,drop=FALSE])}else{0})
            dP_domega <- control$omegajacobian(omegaorig)*dP_domega # this puts the derivative back on the correct scale dL/dpsi = dL/dtheta * dtheta/dpsi, e.g. psi=log(theta)
     
            grad <- c(dP_dbeta,dP_domega)
        }
        else{ #censoringtype=="interval" 
            dP_dbeta <- (if(Utest){colSums(X[notcensored,,drop=FALSE] - dJ_dbeta1[notcensored,,drop=FALSE])}else{0}) + 
                        (if(Rtest){colSums(-dJ_dbeta1[rightcensored,,drop=FALSE])}else{0}) + 
                        (if(Ltest){colSums((S1[leftcensored]/(1-S1[leftcensored]))*dJ_dbeta1[leftcensored,,drop=FALSE])}else{0}) + 
                        (if(Itest){colSums((1/(S1[intervalcensored]-S2[intervalcensored]))*(dJ_dbeta2[intervalcensored,]*S2[intervalcensored]-dJ_dbeta1[intervalcensored,,drop=FALSE]*S1[intervalcensored]))}else{0})
            dh_domega <- matrix(haz$gradh(surv[,"time1"]),ncol=length(omega))
            dP_domega <-    (if(Utest){colSums(dh_domega[notcensored,,drop=FALSE] / h[notcensored] - dJ_domega1[notcensored,,drop=FALSE])}else{0}) + 
                            (if(Rtest){colSums(-dJ_domega1[rightcensored,,drop=FALSE])}else{0}) + 
                            (if(Ltest){colSums((S1[leftcensored]/(1-S1[leftcensored]))*dJ_domega1[leftcensored,,drop=FALSE])}else{0}) + 
                            (if(Itest){colSums((1/(S1[intervalcensored]-S2[intervalcensored]))*(dJ_domega2[intervalcensored,]*S2[intervalcensored]-dJ_domega1[intervalcensored,,drop=FALSE]*S1[intervalcensored]))}else{0})
            dP_domega <- control$omegajacobian(omegaorig)*dP_domega # this puts the derivative back on the correct scale dL/dpsi = dL/dtheta * dtheta/dpsi, e.g. psi=log(theta)            
            
            grad <- c(dP_dbeta,dP_domega)
        }
    
        return(-grad)
    }

}

