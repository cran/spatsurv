##' survspat function
##'
##' A function to run a Bayesian analysis on right censored survial data assuming a proportional hazards model with
##' baseline hazard derived from the exponential model.
##'
##' @param formula see ?flexsurvreg 
##' @param data a SpatialPointsDataFrame object
##' @param dist choice of distribution function for baseline hazard options are: "exp"
##' @param cov.model an object of class covmodel, see ?covmodel
##' @param mcmc.control mcmc control parameters, see ?mcmcpars
##' @param priors an object of class Priors, see ?mcmcPriors
##' @param control additional control parameters, see ?inference.control
##' @return the mcmc output
##' @seealso \link{inference.control}
##' @export


survspat <- function(   formula,
                        data,
                        dist,
                        cov.model,
                        mcmc.control,
                        priors,
                        control=inference.control(gridded=FALSE)){
    
    # initial checks
    if(!inherits(data,"SpatialPointsDataFrame")){
        stop("'data' must be of class 'SpatialPointsDataFrame'.")
    }                    
    responsename <- as.character(formula[[2]])
    survivaldata <- data@data[[responsename]]
    checkSurvivalData(survivaldata)                        

    # okay, start the MCMC!
    start <- Sys.time()

    if(!inherits(data,"SpatialPointsDataFrame")){
        stop("data must be an object of class SpatialPointsDataFrame")
    }
    
    coords <- coordinates(data)

    control$dist <- dist
    funtxt <- ""    
    if(control$gridded){
        funtxt <- "_gridded"
    }
    
    gridobj <- NULL
   
    if(control$gridded){
        gridobj <- FFTgrid(spatialdata=data,cellwidth=control$cellwidth,ext=control$ext)
    	del1 <- gridobj$del1
    	del2 <- gridobj$del2
    	Mext <- gridobj$Mext
    	Next <- gridobj$Next
    	mcens <- gridobj$mcens
    	ncens <- gridobj$ncens    	
    	## COMPUTE GRID DISTANCES ##
    	x <- gridobj$mcens
        y <- gridobj$ncens    
        xidx <- rep(1:Mext,Next)
        yidx <- rep(1:Next,each=Mext)
        dxidx <- pmin(abs(xidx-xidx[1]),Mext-abs(xidx-xidx[1]))
        dyidx <- pmin(abs(yidx-yidx[1]),Next-abs(yidx-yidx[1]))
        u <- sqrt(((x[2]-x[1])*dxidx)^2+((y[2]-y[1])*dyidx)^2)
        
        spix <- grid2spix(xgrid=mcens,ygrid=mcens)
        
        control$fftgrid <- gridobj
        control$idx <- over(data,geometry(spix))
        control$Mext <- Mext
        control$Next <- Next
        control$uqidx <- unique(control$idx)
        
        cat("Output grid size: ",Mext/control$ext," x ",Next/control$ext,"\n")
    }
    else{
        u <- as.vector(as.matrix(dist(coords)))
    }

    DATA <- data    
    
    data <- data@data  
                        
    ##########
    # This chunk of code borrowed from flexsurvreg    
    ##########  
                      
    call <- match.call()
    indx <- match(c("formula", "data"), names(call), nomatch = 0)
    if (indx[1] == 0){ 
        stop("A \"formula\" argument is required")
    }
    temp <- call[c(1, indx)]
    temp[[1]] <- as.name("model.frame")
    m <- eval(temp, parent.frame())
    Y <- model.extract(m, "response")
    if (!inherits(Y, "Surv")){ 
        stop("Response must be a survival object")
    }
    if (!(attr(Y, "type") %in% c("right", "counting"))){ 
        stop("Survival object type \"", attr(Y, "type"), "\"", 
            " not supported")
    }
    if (attr(Y, "type") == "counting"){ 
        Y <- cbind(Y, time = Y[, "stop"] - Y[, "start"])
    }
    else{ 
        Y <- cbind(Y, start = 0, stop = Y[, "time"])
    }
    Terms <- attr(m, "terms")
    X <- model.matrix(Terms, m)
    dat <- list(Y = Y, X = X[, -1, drop = FALSE], Xraw = m[, 
        -1, drop = FALSE])
    X <- dat$X
    
    ##########
    # End of borrowed code    
    ##########                                       
                        
    mcmcloop <- mcmcLoop(N=mcmc.control$nits,burnin=mcmc.control$burn,thin=mcmc.control$thin,progressor=mcmcProgressTextBar)                            
              
                  
                        
    tm <- Y[,1]
    delta <- Y[,2]
    
    cat("\n","Getting initial estimates of model parameters using flexsurvreg","\n")
    mlmod <- suppressWarnings(flexsurvreg(formula,data=data,dist=dist))
    cat("Done.\n")
    estim <- mlmod$opt$par
    print(mlmod)
    #browser() 
    
    cat("Calibrating MCMC algorithm and finding initial values ...\n")
    
    
    
    trns <- getomegatrans(dist)  
    
    if(dist=="exp"){    
        betahat <- estim[2:length(estim)]
        omegahat <- do.call(paste("transformestimates.",dist,sep=""),args=list(x=exp(estim[1]))) 
        omegahat <- log(omegahat)
        omegatrans <- trns$trans
        omegaitrans <- trns$itrans
    }
    else if(dist=="weibull"){    
        betahat <- estim[3:length(estim)]
        omegahat <- do.call(paste("transformestimates.",dist,sep=""),args=list(x=exp(estim[1:2])))
        omegahat <- log(omegahat)
        omegatrans <- trns$trans
        omegaitrans <- trns$itrans
    }
    else{
        stop("Unknown dist, must be one of 'exp' or 'weibull'")    
    }
    
    control$omegatrans <- omegatrans
    control$omegaitrans <- omegaitrans
    control$omegajacobian <- trns$jacobian # used in computing the derivative of the log posterior with respect to the transformed omega (since it is easier to compute with respect to omega) 
    control$omegahessian <- trns$hessian    
    
    control$sigmaidx <- match("sigma",cov.model$parnames)
    if(is.na(control$sigmaidx)){
        stop("At least one of the parameters must be the variance of Y, it should be named sigma")
    }
    
    Yhat <- estimateY(  X=X,
                        betahat=betahat,
                        omegahat=omegahat,
                        surv=survivaldata,
                        control=control)
                        
    calibrate <- get(paste("proposalVariance",funtxt,sep=""))    
       
    other <- calibrate( X=X,
                        surv=survivaldata,
                        betahat=betahat,
                        omegahat=omegahat,
                        Yhat=Yhat,
                        priors=priors,
                        cov.model=cov.model,
                        u=u,
                        control=control) 
    
    #gammahat <- other$gammahat
    etahat <- other$etahat                                                                        
    SIGMA <- other$sigma 

    beta <- betahat
    omega <- omegahat
    eta <- etahat 
 
    gamma <- rep(0,nrow(X))
    if(control$gridded){
        gamma <- matrix(0,control$Mext,control$Next)
    }
        
    lenbeta <- length(beta)
    lenomega <- length(omega)
    leneta <- length(eta)
    lengamma <- length(gamma)
    
    
    
    npars <- lenbeta + lenomega + leneta + lengamma
    
    SIGMA[1:(lenbeta+lenomega),1:(lenbeta+lenomega)] <- (1.65^2/((lenbeta+lenomega)^(1/3)))*SIGMA[1:(lenbeta+lenomega),1:(lenbeta+lenomega)]
    SIGMA[(lenbeta+lenomega+1):(lenbeta+lenomega+leneta),(lenbeta+lenomega+1):(lenbeta+lenomega+leneta)] <- 0.4*(2.38^2/leneta)* SIGMA[(lenbeta+lenomega+1):(lenbeta+lenomega+leneta),(lenbeta+lenomega+1):(lenbeta+lenomega+leneta)]
    SIGMA[(lenbeta+lenomega+leneta+1):(lenbeta+lenomega+leneta+lengamma),(lenbeta+lenomega+leneta+1):(lenbeta+lenomega+leneta+lengamma)] <- (1.65^2/(lengamma^(1/3)))*SIGMA[(lenbeta+lenomega+leneta+1):(lenbeta+lenomega+leneta+lengamma),(lenbeta+lenomega+leneta+1):(lenbeta+lenomega+leneta+lengamma)]   
    
    if(control$gridded){
        matidx <- matrix(0,control$Mext,control$Next)
        matidx[1:(control$Mext/control$ext),1:(control$Next/control$ext)] <- 1
        matidx <- as.logical(matidx) # used to select which Y's to save 
    }
  
    
    diagidx <- 1:npars
    diagidx <- matrix(diagidx,nrow=npars,ncol=2)
    SIGMApars <- as.matrix(SIGMA[1:(lenbeta+lenomega+leneta),1:(lenbeta+lenomega+leneta)])
    SIGMAparsINV <- solve(SIGMApars)
    cholSIGMApars <- t(chol(SIGMApars))  
    #browser()  
    SIGMAgamma <- SIGMA[diagidx][(lenbeta+lenomega+leneta+1):npars]
    SIGMAgammaINV <- 1/SIGMAgamma
    cholSIGMAgamma <- sqrt(SIGMAgamma)
    
     
    
    cat("Running MCMC ...\n")
    
    h <- 1
    
    
    
    LOGPOST <- get(paste("logPosterior",funtxt,sep=""))
    
    oldlogpost <- LOGPOST(  surv=survivaldata,
                            X=X,
                            beta=beta,
                            omega=omega,
                            eta=eta,
                            gamma=gamma,
                            priors=priors,
                            cov.model=cov.model,
                            u=u,
                            control=control,
                            gradient=TRUE)
                            
                          
                                                        
    
    betasamp <- c()
    omegasamp <- c()
    etasamp <- c()
    Ysamp <- c()
    
    gamma <- c(gamma) # turn gamma into a vector 
    
    tarrec <- oldlogpost$logpost
    
    print(SIGMA[1:8,1:8])
    
    
    while(nextStep(mcmcloop)){

        stuffpars <- c(beta,omega,eta)
        propmeanpars <- stuffpars + (h/2)*SIGMApars%*%oldlogpost$grad[1:(lenbeta+lenomega+leneta)]
        newstuffpars <- propmeanpars + sqrt(h)*cholSIGMApars%*%rnorm(lenbeta+lenomega+leneta)
        
 
        propmeangamma <- gamma + (h/2)*SIGMAgamma*oldlogpost$grad[(lenbeta+lenomega+leneta+1):npars]
        newstuffgamma <- propmeangamma + sqrt(h)*cholSIGMAgamma*rnorm(lengamma)
        ngam <- newstuffgamma
        if(control$gridded){
            ngam <- matrix(ngam,control$Mext,control$Next)
        }
                                       
        newlogpost <- LOGPOST(  surv=survivaldata,
                                X=X,
                                beta=newstuffpars[1:lenbeta],
                                omega=newstuffpars[(lenbeta+1):(lenbeta+lenomega)],
                                eta=newstuffpars[(lenbeta+lenomega+1):(lenbeta+lenomega+leneta)],
                                gamma=ngam,
                                priors=priors,
                                cov.model=cov.model,
                                u=u,
                                control=control,
                                gradient=TRUE)

        revmeanpars <- newstuffpars + (h/2)*SIGMApars%*%newlogpost$grad[1:(lenbeta+lenomega+leneta)]
        revmeangamma <- newstuffgamma + (h/2)*SIGMAgamma*newlogpost$grad[(lenbeta+lenomega+leneta+1):npars]       

        revdiffpars <- as.matrix(stuffpars-revmeanpars)
        forwdiffpars <- as.matrix(newstuffpars-propmeanpars)
        revdiffgamma <- as.matrix(gamma-revmeangamma)
        forwdiffgamma <- as.matrix(newstuffgamma-propmeangamma)


        logfrac <- newlogpost$logpost - oldlogpost$logpost - 
                            (0.5/h)*t(revdiffpars)%*%SIGMAparsINV%*%revdiffpars + 
                            (0.5/h)*t(forwdiffpars)%*%SIGMAparsINV%*%forwdiffpars -
                            (0.5/h)*sum(revdiffgamma*SIGMAgammaINV*revdiffgamma) + 
                            (0.5/h)*sum(forwdiffgamma*SIGMAgammaINV*forwdiffgamma)
        
        ac <- min(1,exp(as.numeric(logfrac)))
        
        if(ac>runif(1)){
            beta <- newstuffpars[1:lenbeta]
            omega <- newstuffpars[(lenbeta+1):(lenbeta+lenomega)]
            eta <- newstuffpars[(lenbeta+lenomega+1):(lenbeta+lenomega+leneta)]
            gamma <- newstuffgamma

            oldlogpost <- newlogpost
        }
        
        h <- exp(log(h) + (1/(iteration(mcmcloop)^0.5))*(ac-0.574))
        if(iteration(mcmcloop)%%100==0){
            cat("\n","h =",h,"\n")
        }
        
        
        if(is.retain(mcmcloop)){
            betasamp <- rbind(betasamp,as.vector(beta))
            omegasamp <- rbind(omegasamp,as.vector(omega))
            etasamp <- rbind(etasamp,as.vector(eta))
            if(control$gridded){
                Ysamp <- rbind(Ysamp,as.vector(oldlogpost$Y[matidx]))
            }
            else{
                Ysamp <- rbind(Ysamp,as.vector(oldlogpost$Y))
            }
            tarrec <- c(tarrec,oldlogpost$logpost)
        }
    }
    
    #par(mfrow=c(2,1))
    #matplot(betasamp,type="s")
    #matplot(exp(omegasamp),type="s")

    retlist <- list()
    retlist$formula <- formula
    retlist$data <- DATA
    retlist$dist <- dist
    retlist$cov.model <- cov.model
    retlist$mcmc.control <- mcmc.control
    retlist$priors <- priors
    retlist$control <- control
    
    retlist$terms <- Terms
    retlist$mlmod <- mlmod
    
    ####
    #   Back transform for output
    ####
    
    omegasamp <- omegaitrans(omegasamp)
    itrans <- get(paste("invtransformestimates.",dist,sep=""))
    if(ncol(omegasamp)==1){    
        omegasamp <- matrix(apply(omegasamp,1,itrans))
    }
    else{
        omegasamp <- t(apply(omegasamp,1,itrans))
    }    
    omegasamp <- labelomegamatrix(m=omegasamp,dist=dist)
    
    etasamp <- sapply(1:length(cov.model$itrans),function(i){cov.model$itrans[[i]](etasamp[,i])})
    colnames(etasamp) <- cov.model$parnames     
    
    ####    
    
    colnames(betasamp) <- attr(Terms,"term.labels")
    retlist$betasamp <- betasamp
    retlist$omegasamp <- omegasamp
    retlist$etasamp <- etasamp
    retlist$Ysamp <- Ysamp
    
    retlist$survivaldata <- survivaldata
    
    retlist$gridded <- control$gridded
    if(control$gridded){
        retlist$M <- Mext/control$ext
        retlist$N <- Next/control$ext
        retlist$xvals <- mcens[1:retlist$M]
        retlist$yvals <- ncens[1:retlist$N]
    }
    
    retlist$tarrec <- tarrec
    retlist$lasth <- h
    
    retlist$omegatrans <- omegatrans
    retlist$omegaitrans <- omegaitrans
    
    retlist$time.taken <- Sys.time() - start
    
    cat("Time taken:",retlist$time.taken,"\n")
    
    class(retlist) <- c("list","mcmcspatsurv")

    return(retlist)
}


