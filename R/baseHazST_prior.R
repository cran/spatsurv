## estimation of loglikelihood for prior ST
# gamma vector of length n + 1, n: number of time intervals. gamma_0 is from N(0, 1)
 
##' indepGaussianpriorST function
##'
##' A function to 
##'
##' @param beta X 
##' @param omega X 
##' @param eta X 
##' @param priors X 
##' @return ...
##' @export

indepGaussianpriorST <- function (beta = NULL, omega = NULL, eta = NULL, priors) {

    lp <- 0
    if (!is.null(priors$betaprior)) {
        lp <- lp + sum(dnorm(beta, priors$betaprior$mean, priors$betaprior$sd,log = TRUE))
    }

    if (!is.null(priors$omegaprior)) {
        lp <- lp + priors$omegaprior$eval(omega)
    }

    if (!is.null(priors$etaprior)) {
        lp <- lp + sum(dnorm(eta, priors$etaprior$mean, priors$etaprior$sd, log = TRUE))
    }

    return(lp)
}




##' derivindepGaussianpriorST function
##'
##' A function to 
##'
##' @param beta X 
##' @param omega X 
##' @param eta X 
##' @param priors X 
##' @return ...
##' @export

derivindepGaussianpriorST <- function (beta = NULL, omega = NULL, eta = NULL, priors){
    if (is.null(beta)){
        deriv1beta <- c()
        deriv2beta <- c()
    }
    else{
        deriv1beta <- (-1/priors$betaprior$sd^2) * (beta - priors$betaprior$mean)
        sdbeta <- priors$betaprior$sd
        if (length(priors$betaprior$sd) < length(beta)) {
            sdbeta <- rep(priors$betaprior$sd, length(beta))
        }
        deriv2beta <- -1/sdbeta^2
    }

    if(is.null(omega)){
        deriv1omega <- c()
        deriv2omega <- c()
    }
    else{
        omegaderiv <- priors$omegaprior$deriv(omega)
        deriv1omega <- omegaderiv$deriv1
        deriv2omega <- omegaderiv$deriv2
    }

    if (is.null(eta)){
        deriv1eta <- c()
        deriv2eta <- c()
    }
    else {
        deriv1eta <- (-1/priors$etaprior$sd^2) * (eta - priors$etaprior$mean)
        sdeta <- priors$etaprior$sd
        if (length(priors$etaprior$sd) < length(eta)) {
            sdeta <- rep(priors$etaprior$sd, length(eta))
        }
        deriv2eta <- -1/sdeta^2
    } 

    deriv1 <- c(deriv1beta, deriv1omega, deriv1eta)
    deriv2 <- matrix(0,length(deriv1),length(deriv1))
    diag(deriv2)[1:length(beta)] <- deriv2beta
    deriv2[(length(beta)+1):(length(beta)+length(omega)),(length(beta)+1):(length(beta)+length(omega))] <- deriv2omega
    diag(deriv2)[(length(beta)+length(omega)+1):(length(beta)+length(omega)+length(eta))] <- deriv2eta

    return(list(deriv1 = deriv1, deriv2 = deriv2))
}


## omegaprior: prior for parameter for baseline hazard function f.
## gprior: prior for parameters of function g.
# modify omega prior ST
# parameters for function f in baseline hazard (same as omega in previous setting)
# parameters for function g in baseline hazard (tau, gamma)
# mean of each parameter should have the same length as number of parameters
## define prior for function g (ie. mean for tau and gamma, sd for tau and gamma)


##' omegapriorGaussST function
##'
##' A function to 
##'
##' @param basehaz X 
##' @param fmean X 
##' @param fsd X 
##' @param taumean X 
##' @param tausd X 
##' @param thetamean X 
##' @param thetasd X 
##' @return ...
##' @export

omegapriorGaussST <- function(basehaz, fmean, fsd, taumean, tausd, thetamean, thetasd){
  
    npars <- distinfo(basehaz)()$npars 
    nfpars <- distinfo(basehaz)()$nfpars
    retlist <- list()
    #retlist$mean <- c(fmean, taumean, thetamean)
    retlist$eval <- function(pars){
        f.par <- pars[1 : nfpars]
        tau <- pars[nfpars + 1]
        theta <- pars[nfpars + 2]
        gamma <- pars[(nfpars+3):npars]
        a <- basehaz$a(exp(theta))

        lp <- sum(dnorm(f.par, fmean, fsd, log = TRUE)) + 
                sum(dnorm(tau, taumean, tausd, log = TRUE)) + 
                    sum(dnorm(theta, thetamean, thetasd, log = TRUE)) + 
                        (-0.5)*gamma[1]^2

        #for(i in 2 : length(gamma)){
            lapply(2 : length(gamma),function(i){lp <<- lp - (gamma[i]- a[i - 1] * gamma[i-1])/(2* (1 - a[i - 1]^2))})
        #}

        return(lp)
    }
  
    retlist$deriv <-function(pars){
        pars.f <- pars[1 : nfpars]
        tau <- pars[nfpars + 1]
        theta <- pars[nfpars + 2]
        gamma <- pars[(nfpars+3):npars]
        a <- basehaz$a(exp(theta))

        deriv1 <- c((-1/fsd^2) * (pars.f - fmean), (-1/tausd^2) * (tau - taumean), (-1/thetasd^2) * (theta - thetamean), -gamma[1])
     
        derivgam <- c()
        derivgam2 <- matrix(0, ncol = length(gamma), nrow = length(gamma))
        derivgam2[1,1] <- -1
        
        #for(i in 2 : length(gamma)){
        lapply(2 : length(gamma),function(i){
            derivgam[i - 1] <<- -(gamma[i] - a[i - 1] * gamma[i-1])/(1 - a[i - 1]^2)
            derivgam2[i,i] <<- -1 / (1 - a[i-1]^2)
            derivgam2[i,i-1] <<- a[i-1] / (1 - a[i-1]^2)
            derivgam2[i-1,i] <<- a[i-1] / (1 - a[i-1]^2)
        }
        )
        #}
        deriv1 <- c(deriv1, derivgam)

        deriv2 <- matrix(0,npars,npars)#as.matrix(Matrix::bdiag(-1/fsd^2, -1/tausd^2, -1/thetasd^2, derivgam2))
        diag(deriv2)[1:nfpars] <- -1/fsd^2
        deriv2[nfpars+1,nfpars+1] <- -1/tausd^2
        deriv2[nfpars+2,nfpars+2] <- -1/thetasd^2
        deriv2[(nfpars+3):npars,(nfpars+3):npars] <- derivgam2

        return(list(deriv1=deriv1,deriv2=deriv2))
    }
  
    #retlist$omegaprior <- omegapriorGauss(mean = fmean, sd = fsd)
    #retlist$gprior <- gpriorGaussianST(mean = c(taumean, thetamean), sd = c(tausd, thetasd))

    class(retlist) <- "omegapriorGauss"
    return(retlist)
}



