##' transformestimates.exp function
##'
##' A function to transform estimates of the parameters of the exponential baseline hazard function, so they are commensurate with R's inbuilt density functions.
##'
##' @param x a vector of paramters 
##' @return the transformed parameters. For the exponential model this is just the identity.
##' @export
transformestimates.exp <- function(x){
    return(x) # note this is logged later for use in the MCMC
}



##' invtransformestimates.exp function
##'
##' A function to transform estimates of the parameters of the exponential baseline hazard function, so they are commensurate with R's inbuilt density functions.
##'
##' @param x a vector of paramters 
##' @return the transformed parameters. For the exponential model this is just the identity.
##' @export
invtransformestimates.exp <- function(x){
    return(x)
}




#################################################################################
# exponential survival model
#################################################################################

##' basehazard.exp function
##'
##' A function to 
##'
##' @param pars X 
##' @return ...
##' @export

basehazard.exp <- function(pars){
    fun <- function(t){
        return(rep(pars,length(t))) # in this case pars is a 1-vector, the rate    
    }
    return(fun)
}



##' gradbasehazard.exp function
##'
##' A function to 
##'
##' @param pars X 
##' @return ...
##' @export

gradbasehazard.exp <- function(pars){
    fun <- function(t){
        return(rep(1,length(t))) # in this case pars is a 1-vector, the rate
    }
    return(fun)
}


##' hessbasehazard.exp function
##'
##' A function to 
##'
##' @param pars X 
##' @return ...
##' @export

hessbasehazard.exp <- function(pars){
    fun <- function(t){
        return(as.list(rep(0,length(t)))) # in this case pars is a 1-vector, the rate
    }
    return(fun)
}



##' cumbasehazard.exp function
##'
##' A function to 
##'
##' @param pars X 
##' @return ...
##' @export

cumbasehazard.exp <- function(pars){
    fun <- function(t){
        return(pars*t) # in this case pars is a 1-vector, the rate
    }
    return(fun)  
}



##' gradcumbasehazard.exp function
##'
##' A function to 
##'
##' @param pars X 
##' @return ...
##' @export

gradcumbasehazard.exp <- function(pars){
    fun <- function(t){
        return(t) # in this case pars is a 1-vector, the rate
    }
    return(fun)    
}



##' hesscumbasehazard.exp function
##'
##' A function to 
##'
##' @param pars X 
##' @return ...
##' @export

hesscumbasehazard.exp <- function(pars){
    fun <- function(t){
        return(as.list(rep(0,length(t)))) # in this case pars is a 1-vector, the rate
    }
    return(fun)
}



#################################################################################



##' densityquantile.exp function
##'
##' A function to 
##'
##' @param pars X 
##' @param other X
##' @return ...
##' @export

densityquantile.exp <- function(pars,other){
    fun <- function(probs){
        return(-log(1-probs)/(pars*other$expXbetaplusY))
    }
    return(fun)    
}