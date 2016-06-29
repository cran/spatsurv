 
## Sample gamma from prior distribution
# pars <- c(pars.f, pars.g, theta, gamma)

# h0 = f * exp(g)
# let f be a function in survspat. eg. exponentialHaz
# breaks: a vector of breaking points for intervals, where starting and end points are min(t0) and max(t + t0)
# t: survival times (duration)
# t_0: real entering time for individuals
# gamma to be simulated from MCMC; matrix with rows wrt different time interval index


##' baseHazST function
##'
##' A function to 
##'
##' @param bh1 X 
##' @param survobj X 
##' @param t0 X 
##' @param nbreaks X 
##' @param breakmethod X 
##' @param MLinits X 
##' @return ...
##' @export

baseHazST <- function(bh1 = NULL, survobj, t0, nbreaks = 5, breakmethod = "quantile", MLinits = NULL){

    Sys.sleep(0.2)
    warning("!!!UNDER TESTING, YOU HAVE BEEN WARNED !!!",immediate.=TRUE)
    warning("!!!UNDER TESTING, YOU HAVE BEEN WARNED !!!")
    Sys.sleep(0.2)

    
    if(any(attr(survobj, "type") == c("left", "right"))){
        times <- survobj[,"time"]
        t_plus_t0 <- as.numeric(times + t0)
    }
    else{
        stop("GPHaz only supports right or left censored data at present.")
        #times <- c(survobj[,"time1"],survobj[,"time2"])
        #t_plus_t_0 <- times + rep(t_0,2) # note this is not currently going to work with interval censored data, let's worry about that later.
    }
  
    if(breakmethod == "quantile"){
        brks <- quantile(c(t0,t_plus_t0), prob = seq(0, 1, length.out = nbreaks))
        idx0 <- as.numeric(cut(t0, breaks= brks, include.lowest = TRUE))
        idx <- as.numeric(cut(t_plus_t0, breaks= brks, include.lowest = TRUE))
    }
    else if(breakmethod == "equal"){
        brks <- seq(min(t0), max(t_plus_t0), length.out = nbreaks)
        idx0 <- as.numeric(cut(t0, breaks = brks, include.lowest = TRUE))
        idx <- as.numeric(cut(t_plus_t0, breaks = brks, include.lowest = TRUE))
    }
    else{
        stop("Undefined break method")
    }

    tdiff <- diff(midpts(brks)) #c(Inf, diff(midpts(brks)))
  
    npf <- bh1$distinfo()$npars
    npg <- 2 + nbreaks - 1 # nbreaks - 1 = number of gammas
    np <- npf + npg
  
    getg <- function(idx,pars){ # note that "idx" will be the same length as "t"
        #browser()
        tau <- pars[1]

        gamma <- pars[3 : length(pars)]
        ng <- length(gamma)
        retlist <- list()
      
        retlist$g <- function(t){ # note we'll be passing t+t_0 to this function
        return(tau * gamma[idx] - tau^2/2)
        }
    
        retlist$dg_dtau <- function(t){
            return(gamma[idx] - tau)
        }
    
        retlist$dg_dgamma <- function(t){
            mm <- matrix(0,length(t), ng)
            sapply(1:length(t),function(i){mm[i,idx[i]] <<- tau})
            return(mm)
            #return(matrix(tau, length(t), ng))
        }
    
        retlist$d2g_dtau_dgamma <- function(t){
            return(lapply(1 : length(t), function(x){vv<-rep(0,ng);vv[idx[x]] <- 1;return(vv)}))
            #return(lapply(1 : length(t), function(x){return(rep(1, ng))}))
        }
    
        retlist$d2g_dtau2 <- function(t){
            return(rep(-1,length(t))) # each observation belongs to exactly one interval #lapply(1 : length(t), function(x){return(rep(-1, ng))}))
        }
    
        retlist$d2g_dgamma2 <- function(t){
            return(lapply(1 : length(t), function(x){return(diag(0, ng))}))
        }
        return(retlist)
    }
    
    flist <- list()

    
    flist$distinfo <- function(){
      
        retlist <- list()

        retlist$npars <- np
        retlist$nfpars <- npf
        retlist$parnames <- c(bh1$distinfo()$parnames, "tau", "theta", paste("gamma", 1 : (npg - 2), sep = ""))
        
        retlist$trans <- function(x){
            ans <- c(bh1$distinfo()$trans(x[1 : npf]), log(x[npf + 1]), log(x[npf + 2]), x[(npf + 3) : np])
            return(ans)
        }

        retlist$itrans <- function(x){
            ans <- c(bh1$distinfo()$itrans(x[1 : npf]), exp(x[npf + 1]), exp(x[npf + 2]), x[(npf + 3) : np])
            return(ans)
        }

        retlist$jacobian <- function(x){
            ans <- c(bh1$distinfo()$jacobian(x[1 : npf]), exp(x[npf + 1]), exp(x[npf + 2]), rep(1,(npg - 2)))
            return(ans)
        }
        retlist$hessian <- c(bh1$distinfo()$hessian, exp, exp, lapply(1 : (npg - 2), function(x){function(x){return(0)}}))


        if(is.null(MLinits)){
            retlist$MLinits <- rep(0, np)
            retlist$MLinits[npf+2] <- -10
        }
        else{
            retlist$MLinits <- MLinits
        }

        retlist$MLmethod <- "optifix"
        retlist$MLfixpars <- c(rep(FALSE,npf+1),rep(TRUE,npg-1))

        return(retlist)

    }


    flist$a <- function(pars){
        return(exp(-pars * tdiff))
    }

   
    flist$basehazard <- function(pars){
        f <- basehazard(bh1)(pars[1 : npf])
        gfuns <- getg(idx, pars[(npf + 1) : np])
        fun <- function(t){     
            return(f(t) * exp(gfuns$g(t+t0))) # as far as g goes, t is a dummy
        }
        return(fun)
    }
    
    flist$gradbasehazard <- function(pars){
        f <- basehazard(bh1)(pars[1:npf])
        fdash <- gradbasehazard(bh1)(pars[1:npf])
        gfuns <- getg(idx, pars[(npf+1):np])
        fun <- function(t){         
            return(cbind(fdash(t), f(t) * cbind(gfuns$dg_dtau(t+t0),0, gfuns$dg_dgamma(t+t0))) * exp(gfuns$g(t+t0))) # note dh/dtheta = 0 ... the derivative w.r.to theta plays a part in the prior
        }
        return(fun)
    }
    
    flist$hessbasehazard <- function(pars){
        f <- basehazard(bh1)(pars[1 : npf])
        fdash <- gradbasehazard(bh1)(pars[1 : npf])
        fhess <- hessbasehazard(bh1)(pars[1 : npf])
        gfuns <- getg(idx, pars[(npf + 1) : np])
        
        fun <- function(t){
            
            g <- gfuns$g(t+t0)
            gradtau <- gfuns$dg_dtau(t+t0)
            gradgam <- gfuns$dg_dgamma(t+t0)
            grad2tau <- gfuns$d2g_dtau2(t+t0)
            grad2gam <- gfuns$d2g_dgamma2(t+t0)
            grad2taugam <- gfuns$d2g_dtau_dgamma(t+t0)
            
            f.base <- f(t)
            f.grad <- matrix(fdash(t),ncol=npf)
            f.hess <- fhess(t)

            expg <- exp(g)
            fexpg <- f.base * expg

            gradtau_gradomega <- lapply(1:length(t),function(i){as.vector(outer(gradtau[i],f.grad[i,]))})
            gradgamma_gradomega <- lapply(1:length(t),function(i){outer(gradgam[i,],f.grad[i,])})
            gradgamma_gradtau <- lapply(1:length(t),function(i){as.vector(outer(gradgam[i,],gradtau[i]))})
            gradgamma_gradgamma <- lapply(1:length(t),function(i){outer(gradgam[i,],gradgam[i,])})

            hess_omega <- mapply("*",f.hess,expg)
            hess_tau <- fexpg*(grad2tau+gradtau^2)
            hess_tau_omega <- mapply("*",gradtau_gradomega,expg)
            hess_gamma <- mapply("+",mapply("*",fexpg,grad2gam),mapply("*",fexpg,gradgamma_gradgamma))
            hess_gamma_omega <- mapply("*",gradgamma_gradomega,expg)
            hess_gamma_tau <- mapply("+",mapply("*",fexpg,grad2taugam),mapply("*",fexpg,gradgamma_gradtau))
            
            funfun <- function(i){
                mat <- matrix(0,np,np)
                mat[1:npf,1:npf] <- hess_omega[[i]]
                mat[1:npf,npf+1] <- mat[npf+1,1:npf] <- hess_tau_omega[[i]]
                mat[npf+1,npf+1] <- hess_tau[[i]]
                mat[(npf+3):np,(npf+3):np] <- hess_gamma[[i]]
                mat[(npf+3):np,1:npf] <- hess_gamma_omega[[i]]
                mat[1:npf,(npf+3):np] <- t(hess_gamma_omega[[i]])
                mat[npf+1,(npf+3):np] <- mat[(npf+3):np,npf+1] <- hess_gamma_tau[[i]]
                return(mat)
            }
            
            hess <- lapply(1:length(t),funfun)
            
            return(hess)
        
        }

        return(fun)
          
    }

    ####################

    idxls <- lapply(1:length(idx),function(i){idx0[i]:idx[i]})
    diffbrks <- diff(brks)
    timeends <- function(i){
        if(length(idxls[[i]])==1){
            tms <- c(0,times[i])
        }
        else if(length(idxls[[i]])==2){
            tms <- c(0,brks[idxls[[i]][1]+1] - t0[i],times[i])
        }
        else{
            tms <- brks[idxls[[i]]][-1] - t0[i]
            tms <- c(0,tms,times[i])
        }
        return(tms)
    }
    tends <- lapply(1:length(idxls),timeends)

    ####################

    flist$cumbasehazard <- function(pars){
        
        F <- cumbasehazard(bh1)(pars[1 : npf])

        tau <- pars[npf+1]
        tau2 <- tau^2
        gamma <- pars[(npf+3):np]

        g <- -tau2/2 + tau*gamma
        expg <- exp(g)
          
        fun <- function(t){
            Fbreaks <- lapply(tends,F)

            deltaF <- lapply(Fbreaks,diff)
            H <- sapply(1:length(idx),function(i){sum(expg[idxls[[i]]]*deltaF[[i]])})

            return(H)
        }
        
        return(fun)
        
    }
    

    flist$gradcumbasehazard <- function(pars){
        
        F <- cumbasehazard(bh1)(pars[1 : npf])
        gradF <- gradcumbasehazard(bh1)(pars[1 : npf])

        tau <- pars[npf+1]
        tau2 <- tau^2
        gamma <- pars[(npf+3):np]

        g <- -tau2/2 + tau*gamma
        expg <- exp(g)

        dexpg_dtau <- (-tau+gamma)*expg
        dexpg_dgamma <- tau*expg
          
        fun <- function(t){
            Fbreaks <- lapply(tends,F)

            gradFwrap <- function(t){
                ans <- matrix(gradF(t),ncol=npf)
                if(any(t==0)){
                    ans[t==0,] <- 0 
                }
                return(ans)
            }
            gradFbreaks <- lapply(1:length(idxls),function(i){gradFwrap(tends[[i]])})

            graddiff <- function(mat){
                m <- matrix(NA,nrow=nrow(mat)-1,ncol=ncol(mat))
                sapply(2:nrow(mat),function(i){m[i-1,] <<- mat[i,] - mat[i-1,]})
                return(m)
            }

            deltaF <- lapply(Fbreaks,diff)
            deltagradF <- lapply(gradFbreaks,graddiff)

            expgcontrib <- lapply(1:length(idx),function(i){expg[idxls[[i]]]})
            n_expgcontrib <- sapply(expgcontrib,length)
            dexpg_dgamma_contrib <- lapply(1:length(idx),function(i){dexpg_dgamma[idxls[[i]]]})

            gradHomega <- matrix(t(sapply(1:length(idx),function(i){colSums(matrix(expgcontrib[[i]]*deltagradF[[i]],nrow=n_expgcontrib[[i]]))})),nrow=length(idx))
            
            gradHtau <- sapply(1:length(idx),function(i){sum(dexpg_dtau[idxls[[i]]]*deltaF[[i]])})
            
            gradHgamma <- t(sapply(1:length(idx),function(i){out <- rep(0,nbreaks-1);out[idxls[[i]]] <- dexpg_dgamma_contrib[[i]]*deltaF[[i]];return(out)}))

            gradH <- cbind(gradHomega,gradHtau,0,gradHgamma)

            return(gradH)

        }
        
        return(fun)
        
    }


    flist$hesscumbasehazard <- function(pars){
        
        F <- cumbasehazard(bh1)(pars[1 : npf])
        gradF <- gradcumbasehazard(bh1)(pars[1 : npf])
        hessF <- hesscumbasehazard(bh1)(pars[1 : npf])

        tau <- pars[npf+1]
        tau2 <- tau^2
        gamma <- pars[(npf+3):np]

        g <- -tau2/2 + tau*gamma
        expg <- exp(g)

        d2expg_dtau2 <- (-1+(-tau+gamma)^2)*expg
        d2expg_dgamma_dtau <- (1+(-tau+gamma)*tau)*expg
        d2expg_domega_dtau <- (-tau+gamma)*expg
        d2expg_domega_dgamma <- tau*expg
        d2expg_dgamma2 <- tau^2*expg

          
        fun <- function(t){

            Fbreaks <- lapply(tends,F)
            gradFwrap <- function(t){
                ans <- matrix(gradF(t),ncol=npf)
                if(any(t==0)){
                    ans[t==0,] <- 0 
                }
                return(ans)
            }
            gradFbreaks <- lapply(1:length(idxls),function(i){gradFwrap(tends[[i]])})

            hessFwrap <- function(t){
                ans <- hessF(t)
                ind <- which(t==0)
                if(length(ind)>0){
                    for(i in ind){
                        ans[[i]] <- matrix(0,npf,npf)
                    }
                }
                return(ans)
            }
            hessFbreaks <- lapply(1:length(idxls),function(i){hessFwrap(tends[[i]])})

            graddiff <- function(mat){
                m <- matrix(NA,nrow=nrow(mat)-1,ncol=ncol(mat))
                sapply(2:nrow(mat),function(i){m[i-1,] <<- mat[i,] - mat[i-1,]})
                return(m)
            }

            hessdiff <- function(matlist){
                l <- list()
                sapply(2:length(matlist),function(i){l[[i-1]] <<- matlist[[i]] - matlist[[i-1]]})
                return(l)
            }

            deltaF <- lapply(Fbreaks,diff)
            deltagradF <- lapply(gradFbreaks,graddiff)
            deltahessF <- lapply(hessFbreaks,hessdiff)

            expgcontrib <- lapply(1:length(idx),function(i){expg[idxls[[i]]]})
            n_expgcontrib <- sapply(expgcontrib,length)
            d2expg_dtau2_contrib <- lapply(1:length(idx),function(i){d2expg_dtau2[idxls[[i]]]})
            d2expg_dgamma_dtau_contrib <- lapply(1:length(idx),function(i){d2expg_dgamma_dtau[idxls[[i]]]})
            d2expg_domega_dtau_contrib <- lapply(1:length(idx),function(i){d2expg_domega_dtau[idxls[[i]]]})
            d2expg_domega_dgamma_contrib <- lapply(1:length(idx),function(i){d2expg_domega_dgamma[idxls[[i]]]})
            d2expg_dgamma2_contrib <- lapply(1:length(idx),function(i){d2expg_dgamma2[idxls[[i]]]})

            hess_omega <- lapply(1:length(idx),function(i){Reduce("+",mapply("*",deltahessF[[i]],expgcontrib[[i]],SIMPLIFY=FALSE))})
            hess_tau <- lapply(1:length(idx),function(i){sum(d2expg_dtau2[idxls[[i]]]*deltaF[[i]])})
            hess_tau_omega <- lapply(1:length(idx),function(i){colSums(matrix(d2expg_domega_dtau_contrib[[i]]*deltagradF[[i]],nrow=n_expgcontrib[i]))})
            hess_gamma <- lapply(1:length(idx),function(i){dg <- rep(0,nbreaks-1);dg[idxls[[i]]] <- d2expg_dgamma2_contrib[[i]]*deltaF[[i]];return(diag(dg))})
            hess_gamma_omega <- lapply(1:length(idx),function(i){mat <- matrix(0,nbreaks-1,npf);mat[idxls[[i]],] <- d2expg_domega_dgamma_contrib[[i]]*deltagradF[[i]];return(mat)})
            hess_gamma_tau <- lapply(1:length(idx),function(i){vc <- rep(0,nbreaks-1);vc[idxls[[i]]] <- d2expg_dgamma_dtau_contrib[[i]]*deltaF[[i]];return(vc)})
            
            funfun <- function(i){
                mat <- matrix(0,np,np)
                mat[1:npf,1:npf] <- hess_omega[[i]]
                mat[1:npf,npf+1] <- mat[npf+1,1:npf] <- hess_tau_omega[[i]]
                mat[npf+1,npf+1] <- hess_tau[[i]]
                mat[(npf+3):np,(npf+3):np] <- hess_gamma[[i]]
                mat[(npf+3):np,1:npf] <- hess_gamma_omega[[i]]
                mat[1:npf,(npf+3):np] <- t(hess_gamma_omega[[i]])
                mat[npf+1,(npf+3):np] <- mat[(npf+3):np,npf+1] <- hess_gamma_tau[[i]]
                return(mat)
            }
            
            hess <- lapply(1:length(t),funfun)

            return(hess)

        }
        
        return(fun)
        
    }

    
    
    flist$densityquantile <- function(pars,other){
        fun <- function(probs){
            stop("densityquantile not available yet")
          
        return(fun)
        }
    }
    
    
    class(flist) <- c("basehazardspec","list")
    return(flist)
        
      
        
}
      



##' gamma2risk function
##'
##' A function to 
##'
##' @param mod X 
##' @return ...
##' @export
      
gamma2risk <- function(mod){
    idx <- grep("gamma",colnames(mod$omegasamp))
    tau <- mod$omegasamp[,"tau"]
    gamma <- mod$omegasamp[,idx]
    Y <- matrix(-tau^2/2,nrow=nrow(gamma),ncol=ncol(gamma)) + tau*gamma
    return(exp(Y))
}



##' boxplotRisk function
##'
##' A function to 
##'
##' @param g2r X 
##' @return ...
##' @export


boxplotRisk <- function(g2r){
    d <- as.data.frame(cbind(g=as.vector(g2r),t=rep(1:ncol(g2r),each=nrow(g2r))))
    boxplot(d$g~d$t)
}     
      
      
      
      
      