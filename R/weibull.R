##' transformestimates.weibull function
##'
##' A function to transform estimates of the parameters of the weibull baseline hazard function, so they are commensurate with R's inbuilt density functions.
##'
##' @param x a vector of paramters
##' @return the transformed parameters. For the weibull model, this transforms 'shape' 'scale' (see ?dweibull) to 'alpha' and 'lambda' for the MCMC
##' @export
transformestimates.weibull <- function(x){
    a <- x[1] # shape
    b <- x[2] # scale
    alpha <- a
    lambda <- (1/b)^a

    ans <- c(alpha=alpha,lambda=lambda) # note this is logged later for use in the MCMC
    
    return(ans)
}



##' invtransformestimates.weibull function
##'
##' A function to transform estimates of the parameters of the weibull baseline hazard function, so they are commensurate with R's inbuilt density functions.
##'
##' @param x a vector of paramters
##' @return the transformed parameters. For the weibull model, this is the back-transform from 'alpha' and 'lambda' to 'shape' 'scale' (see ?dweibull).
##' @export
invtransformestimates.weibull <- function(x){

    alpha <- x[1]
    lambda <- x[2]
    
    shape <- alpha
    scale <- exp((-1/alpha)*log(lambda))

    ans <- c(shape=shape,scale=scale)    
    
    return(ans)
}



#################################################################################
# Weibull survival model
#################################################################################


##' basehazard.weibull function
##'
##' A function to 
##'
##' @param pars X 
##' @return ...
##' @export

basehazard.weibull <- function(pars){
    fun <- function(t){
        return(pars[2]*pars[1]*t^(pars[1]-1)) # in this case alpha=pars[1], lambda=pars[2] 
    }
    return(fun)  
}



##' gradbasehazard.weibull function
##'
##' A function to 
##'
##' @param pars X 
##' @return ...
##' @export

gradbasehazard.weibull <- function(pars){
    fun <- function(t){
        return(t^(pars[1]-1)*cbind(pars[2]*(1+pars[1]*log(t)),pars[1])) # in this case alpha=pars[1], lambda=pars[2]
    }
    return(fun)
    
}



##' hessbasehazard.weibull function
##'
##' A function to 
##'
##' @param pars X 
##' @return ...
##' @export

hessbasehazard.weibull <- function(pars){
    funfun <- function(t,pars){
        m <- matrix(0,2,2) # note m[2,2]=0 i.e. d2h_0/dlambda^2 = 0
        m[1,2] <- m[2,1] <- t^(pars[1]-1)*(1+pars[1]*log(t))
        m[1,1] <- pars[2]*t^(pars[1]-1)*log(t)*(2+pars[1]*log(t))
        return(as.vector(m)) # in this case alpha=pars[1], lambda=pars[2]
    }
    
    fun <- function(t){
        return(lapply(t,funfun,pars=pars))
    }
    return(fun)
    
}



##' cumbasehazard.weibull function
##'
##' A function to 
##'
##' @param pars X 
##' @return ...
##' @export

cumbasehazard.weibull <- function(pars){
    fun <- function(t){
        return(pars[2]*t^(pars[1])) # in this case alpha=pars[1], lambda=pars[2]
    }
    return(fun) 
}




##' gradcumbasehazard.weibull function
##'
##' A function to 
##'
##' @param pars X 
##' @return ...
##' @export

gradcumbasehazard.weibull <- function(pars){
    fun <- function(t){
        return(t^(pars[1])*cbind(pars[2]*log(t),1)) # in this case alpha=pars[1], lambda=pars[2]
    }
    return(fun)    
}




##' hesscumbasehazard.weibull function
##'
##' A function to 
##'
##' @param pars X 
##' @return ...
##' @export
hesscumbasehazard.weibull <- function(pars){
    funfun <- function(t,pars){
        m <- matrix(0,2,2) # note m[2,2]=0 i.e. d2H_0/dlambda^2 = 0
        other <- log(t)*t^pars[1]
        m[1,2] <- m[2,1] <- other 
        m[1,1] <- pars[2]*other*log(t)
        return(m) # in this case alpha=pars[1], lambda=pars[2]
    }
    
    fun <- function(t){
        return(lapply(t,funfun,pars=pars))
    }
    return(fun)
}



#################################################################################



##' densityquantile.weibull function
##'
##' A function to 
##'
##' @param pars X 
##' @param other X
##' @return ...
##' @export

densityquantile.weibull <- function(pars,other){
    fun <- function(probs){
        return((-log(1-probs)/(pars[2]*other$expXbetaplusY))^(1/pars[1]))
    }
    return(fun)    
}