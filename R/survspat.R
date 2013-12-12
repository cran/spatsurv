##' survspat function
##'
##' A function to run a Bayesian analysis on right censored survial data assuming a proportional hazards model with
##' baseline hazard derived from the exponential model.
##'
##' @param formula see ?flexsurvreg 
##' @param data a SpatialPointsDataFrame object
##' @param dist choice of distribution function for baseline hazard options are: "exp"
##' @param covmodel an object of class covmodel, see ?covmodel
##' @param mcmc.control mcmc control parameters, see ?mcmcpars
##' @param priors an object of class Priors, see ?mcmcPriors
##' @param control additional control parameters, see ?inference.control
##' @return the mcmc output
##' @seealso \link{inference.control}
##' @export


survspat <- function(formula,data,dist,covmodel,mcmc.control,priors,control=inference.control(gridded=FALSE)){

    start <- Sys.time()

    if(!inherits(data,"SpatialPointsDataFrame")){
        stop("data must be an object of class SpatialPointsDataFrame")
    }
    
    coords <- coordinates(data)

    funtxt <- ""    
    if(control$gridded){
        funtxt <- ".gridded"
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
    
    

    if(dist=="exp"){    
        betahat <- estim[2:length(estim)]
        omegahat <- do.call(paste("transformestimates.",dist,sep=""),args=list(x=exp(estim[1]))) 
        omegahat <- log(omegahat)
    }
    else if(dist=="weibull"){    
        betahat <- estim[3:length(estim)]
        omegahat <- do.call(paste("transformestimates.",dist,sep=""),args=list(x=exp(estim[1:2])))
        omegahat <- log(omegahat) 
    }
    else{
        stop("Unknown dist, must be one of 'exp' or 'weibull'")    
    }
    
    Yhat <- do.call(paste("estimateY.",dist,sep=""),args=list(X=X,betahat=betahat,omegahat=omegahat,tm=tm,delta=delta))    
       
    other <- do.call(paste("proposalvariance.",dist,funtxt,sep=""),args=list(   X=X,
                                                                                delta=delta,
                                                                                tm=tm,
                                                                                betahat=betahat,
                                                                                omegahat=omegahat,
                                                                                Yhat=Yhat,
                                                                                priors=priors,
                                                                                covmodel=covmodel,
                                                                                u=u,
                                                                                control=control)) 
    
    gammahat <- other$gammahat
    etahat <- other$etahat                                                                        
    SIGMA <- other$sigma 

    beta <- betahat
    omega <- omegahat
    eta <- etahat 
 
    gamma <- rep(0,length(gammahat))
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
    
        matidx <- (lenbeta+lenomega+leneta+1):npars
        matidx <- matrix(matidx,nrow=length(matidx),ncol=2)
               
        SIGMAINV <- Matrix(0,npars,npars)
        SIGMAINV[1:(lenbeta+lenomega+leneta),1:(lenbeta+lenomega+leneta)] <- solve(as.matrix(SIGMA[1:(lenbeta+lenomega+leneta),1:(lenbeta+lenomega+leneta)]))
        SIGMAINV[matidx] <- 1/SIGMA[matidx]
        
        cholSIGMA <- Matrix(0,npars,npars)
        cholSIGMA[1:(lenbeta+lenomega+leneta),1:(lenbeta+lenomega+leneta)] <- chol(as.matrix(SIGMA[1:(lenbeta+lenomega+leneta),1:(lenbeta+lenomega+leneta)]))
        cholSIGMA[matidx] <- sqrt(SIGMA[matidx])
        
        matidx <- matrix(0,control$Mext,control$Next)
        matidx[1:(control$Mext/control$ext),1:(control$Next/control$ext)] <- 1
        matidx <- as.logical(matidx)
    }
    else{
        browser()
        SIGMAINV <- solve(SIGMA) # SIGMA is sparse, so this is easy to compute    
        cholSIGMA <- Matrix(t(chol(SIGMA)))
    }
    #print(SIGMA[1:(lenbeta+lenomega),1:(lenbeta+lenomega)])    
    
    cat("Running MCMC ...\n")
    
    h <- 1
    
    
    
    LOGPOST <- get(paste("logposterior.",dist,funtxt,sep=""))
    
    oldlogpost <- LOGPOST(  tm=tm,
                            delta=delta,
                            X=X,beta=beta,
                            omega=omega,
                            eta=eta,
                            gamma=gamma,
                            priors=priors,
                            covmodel=covmodel,
                            u=u,
                            control=control)
                                                        
    
    betasamp <- c()
    omegasamp <- c()
    etasamp <- c()
    Ysamp <- c()
    
    tarrec <- oldlogpost$logpost
    
    print(SIGMA[1:8,1:8])
    
    
    while(nextStep(mcmcloop)){
    
        stuff <- c(beta,omega,eta,gamma)
        propmean <- stuff + (h/2)*SIGMA%*%oldlogpost$grad
        newstuff <- propmean + h*cholSIGMA%*%rnorm(npars)
        
        ngam <- newstuff[(lenbeta+lenomega+leneta+1):npars]
        if(control$gridded){
            ngam <- matrix(ngam,control$Mext,control$Next)
        }
        
        newlogpost <- LOGPOST(  tm=tm,
                                delta=delta,
                                X=X,
                                beta=newstuff[1:lenbeta],
                                omega=newstuff[(lenbeta+1):(lenbeta+lenomega)],
                                eta=newstuff[(lenbeta+lenomega+1):(lenbeta+lenomega+leneta)],
                                gamma=ngam,
                                priors=priors,
                                covmodel=covmodel,
                                u=u,
                                control=control)
                                
        revmean <- newstuff +  (h/2)*SIGMA%*%newlogpost$grad  
         
        
        revdiff <- as.matrix(stuff-revmean)
        forwdiff <- as.matrix(newstuff-propmean)
                       
        logfrac <- newlogpost$logpost - oldlogpost$logpost -0.5*t(revdiff)%*%SIGMAINV%*%revdiff + 0.5*t(forwdiff)%*%SIGMAINV%*%forwdiff
        
        ac <- min(1,exp(as.numeric(logfrac)))
        
        if(ac>runif(1)){
            beta <- newstuff[1:lenbeta]
            omega <- newstuff[(lenbeta+1):(lenbeta+lenomega)]
            eta <- newstuff[(lenbeta+lenomega+1):(lenbeta+lenomega+leneta)]
            gamma <- newstuff[(lenbeta+lenomega+leneta+1):npars]
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
    retlist$formula
    retlist$data
    retlist$mcmc.control
    retlist$betapriormean
    retlist$betapriorsd
    retlist$omegapriormean
    retlist$omegapriorsd

    retlist$terms <- Terms
    retlist$mlmod <- mlmod    
    
    retlist$betasamp <- betasamp
    retlist$omegasamp <- omegasamp
    retlist$etasamp <- etasamp
    retlist$Ysamp <- Ysamp
    
    if(control$gridded){
        retlist$M <- Mext/control$ext
        retlist$N <- Next/control$ext
        retlist$xvals <- mcens[1:retlist$M]
        retlist$yvals <- ncens[1:retlist$N]
    }
    
    retlist$tarrec <- tarrec
    retlist$lasth <- h
    
    retlist$time.taken <- Sys.time() - start
    
    cat("Time taken:",retlist$time.taken,"\n")
    
    class(retlist) <- c("list","mcmcspatsurv")

    return(retlist)
}


