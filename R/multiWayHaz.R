##' multiWayHaz function
##'
##' A function to
##'
##' @param bhlist X
##' @param bhtime X
##' @param bhfix X
##' @param MLinits X
##' @return ...
##' @export

multiWayHaz <- function(bhlist,bhtime,bhfix,MLinits=NULL){

    nhaz <- length(bhlist) # first remove unnecessary parameters from
    hesslist <- list()
    for(i in 1:nhaz){
        if(! (identical(bhlist[[i]]$distinfo()$trans,log) & identical(bhlist[[i]]$distinfo()$itrans,exp) & identical(bhlist[[i]]$distinfo()$jacobian,exp))){
            stop("Check transformations in multiWayHaz ... :-(")
        }
        if(!is.null(bhfix[[i]])){
            bhlist[[i]] <- fixParHaz(bhlist[[i]],bhfix[[i]]$idx,bhfix[[i]]$fixval)
        }
        hesslist <- c(hesslist,bhlist[[i]]$distinfo()$hessian)
    }

    pn <- unlist(sapply(bhlist,function(x){x$distinfo()$parnames}))
    np <- length(pn)
    npvec <- unlist(sapply(bhlist,function(x){x$distinfo()$npars}))
    csnpvec <- cumsum(npvec)
    vf <- function(i){
        if(i ==1){
            return(1:csnpvec[i])
        }
        else{
            return((csnpvec[i-1]+1):csnpvec[i])
        }
    }
    paridx <- lapply(1:length(csnpvec),vf) # indices of parameters to feed to each bhlist function

    post_processing_mode <- FALSE

    flist <- list()

    flist$distinfo <- function(){
        retlist <- list()
        retlist$npars <- np
        retlist$parnames <- pn
        retlist$trans <- log
        retlist$itrans <- exp
        retlist$jacobian <- exp
        retlist$hessian <- hesslist
        retlist$MLinits <- MLinits
        return(retlist)
    }

    flist$basehazard <- function(pars){
        fun <- function(t,...){
            #browser()
            ans <- sapply(1:nhaz,function(i){bhlist[[i]]$basehazard(pars[paridx[[i]]])(bhtime[[i]])})
            if(post_processing_mode){
                return(ans)
            }
            return(apply(ans,1,prod))
        }
        return(fun)
    }

    flist$gradbasehazard <- function(pars){
        fun <- function(t,...){
            h <- sapply(1:nhaz,function(i){bhlist[[i]]$basehazard(pars[paridx[[i]]])(bhtime[[i]])})
            #hdash <- matrix(unlist(sapply(1:nhaz,function(i){apply(h[,-i,drop=FALSE],1,prod)*bhlist[[i]]$gradbasehazard(pars[paridx[[i]]])(bhtime[[i]])})),ncol=np)
            prodh <- apply(h,1,prod)
            hdash <- matrix(unlist(sapply(1:nhaz,function(i){(prodh/h[,i])*bhlist[[i]]$gradbasehazard(pars[paridx[[i]]])(bhtime[[i]])})),ncol=np)
            return(hdash)
        }
        return(fun)
    }

    flist$hessbasehazard <- function(pars){
        fun <- function(t,...){
            h <- sapply(1:nhaz,function(i){bhlist[[i]]$basehazard(pars[paridx[[i]]])(bhtime[[i]])})
            hdash <- matrix(unlist(sapply(1:nhaz,function(i){apply(h[,-i,drop=FALSE],1,prod)*bhlist[[i]]$gradbasehazard(pars[paridx[[i]]])(bhtime[[i]])})),ncol=np)
            hess <- lapply(1:nhaz,function(i){bhlist[[i]]$hessbasehazard(pars[paridx[[i]]])(bhtime[[i]])}) # here and above 2 lines, i is bh index
            getMlist <- function(i){ # here i is subject index
                return(lapply(1:nhaz,function(j){hess[[j]][[i]]}))
            }
            hesslist <- list()
            hmult <- matrix(0,nrow(h),ncol(h))
            for(j in 1:ncol(hmult)){
                hmult[,j] <- apply(h[,-j,drop=FALSE],1,prod)
            }
            for(i in 1:length(bhtime[[1]])){ # loop over subjects, yes this is likely slow, but does not matter ...
                m <- blockDiag(getMlist(i))
                m <- rep(hmult[i,],npvec)*m # scale diagonal entries
                for(j in 1:nhaz){
                    for(k in 1:nhaz){
                        if(j == k){
                            next
                        }
                        else{
                            m[paridx[[j]],paridx[[k]]] <- prod(h[i,-c(j,k)]) * outer(hdash[i,paridx[[j]]],hdash[i,paridx[[k]]])
                        }
                    }
                }
                hesslist[[i]] <- m
            }
            return(hesslist)
        }
        return(fun)

    }

    flist$cumbasehazard <- function(pars){
        fun <- function(t,...){
            ans <- sapply(1:nhaz,function(i){bhlist[[i]]$cumbasehazard(pars[paridx[[i]]])(bhtime[[i]])})
            if(post_processing_mode){
                return(ans)
            }
            return(apply(ans,1,prod))
        }
        return(fun)
    }

    flist$gradcumbasehazard <- function(pars){
        fun <- function(t,...){
            H <- sapply(1:nhaz,function(i){bhlist[[i]]$cumbasehazard(pars[paridx[[i]]])(bhtime[[i]])})
            #hdash <- matrix(unlist(sapply(1:nhaz,function(i){apply(H[,-i,drop=FALSE],1,prod)*bhlist[[i]]$gradcumbasehazard(pars[paridx[[i]]])(bhtime[[i]])})),ncol=np)
            prodH <- apply(H,1,prod)
            hdash <- matrix(unlist(sapply(1:nhaz,function(i){(prodH/H[,i])*bhlist[[i]]$gradcumbasehazard(pars[paridx[[i]]])(bhtime[[i]])})),ncol=np)
            return(hdash)
        }
        return(fun)
    }

    flist$hesscumbasehazard <- function(pars){
        fun <- function(t,...){
            h <- sapply(1:nhaz,function(i){bhlist[[i]]$cumbasehazard(pars[paridx[[i]]])(bhtime[[i]])})
            hdash <- matrix(unlist(sapply(1:nhaz,function(i){apply(h[,-i,drop=FALSE],1,prod)*bhlist[[i]]$gradcumbasehazard(pars[paridx[[i]]])(bhtime[[i]])})),ncol=np)
            hess <- lapply(1:nhaz,function(i){bhlist[[i]]$hesscumbasehazard(pars[paridx[[i]]])(bhtime[[i]])}) # here and above 2 lines, i is bh index
            getMlist <- function(i){ # here i is subject index
                return(lapply(1:nhaz,function(j){hess[[j]][[i]]}))
            }
            hesslist <- list()
            hmult <- matrix(0,nrow(h),ncol(h))
            for(j in 1:ncol(hmult)){
                hmult[,j] <- apply(h[,-j,drop=FALSE],1,prod)
            }
            for(i in 1:length(bhtime[[1]])){ # loop over subjects, yes this is likely slow, but does not matter ...
                m <- blockDiag(getMlist(i))
                m <- rep(hmult[i,],npvec)*m # scale diagonal entries
                for(j in 1:nhaz){
                    for(k in 1:nhaz){
                        if(j == k){
                            next
                        }
                        else{
                            m[paridx[[j]],paridx[[k]]] <- prod(h[i,-c(j,k)]) * outer(hdash[i,paridx[[j]]],hdash[i,paridx[[k]]])
                        }
                    }
                }
                hesslist[[i]] <- m
            }
            return(hesslist)
        }
        return(fun)

    }

    flist$densityquantile <- function(pars,other){
        fun <- function(probs,...){
            return((-log(1-probs)/(pars[2]*other$expXbetaplusY))^(1/pars[1]))
        }
        return(fun)
    }


    class(flist) <- c("basehazardspec","list")
    return(flist)
}


##' insert function
##'
##' A function to
##'
##' @param pars X
##' @param idx X
##' @param val X
##' @return ...
##' @export

insert <- function(pars,idx,val){
    n <- length(pars) # note length(pars) is already reduced by 1 here, so true parameter length is n+1
    if(idx==1){
        return(c(val,pars))
    }
    else if(idx==(n+1)){
        return(c(pars,val))
    }
    else{
        return(c(pars[1:(idx-1)],val,pars[idx:n]))
    }
}


##'
##' A function to
##'
##' @param matlist X
##' @return ...
##' @export

blockDiag <- function(matlist){
    d <- sapply(matlist,function(x){if(is.matrix(x)){return(dim(x)[1])}else{return(1)}})
    n <- sum(d)
    csnpvec <- cumsum(d)
    vf <- function(i){
        if(i ==1){
            return(1:csnpvec[i])
        }
        else{
            return((csnpvec[i-1]+1):csnpvec[i])
        }
    }
    paridx <- lapply(1:length(csnpvec),vf)
    m <- matrix(0,n,n)
    lapply(1:length(matlist),function(i){m[paridx[[i]],paridx[[i]]] <<- matlist[[i]]})
    return(m)
}

##' fixParHaz function
##'
##' A function to
##'
##' @param bh X
##' @param idx X
##' @param fixval X
##' @return ...
##' @export

fixParHaz <- function(bh,idx,fixval){ # accepts objects inheriting class "basehazardspec"

    flist <- list()

    if(! (identical(bh$distinfo()$trans,log) & identical(bh$distinfo()$itrans,exp) & identical(bh$distinfo()$jacobian,exp))){
        stop("Check transformations in fixParHaz ... :-(")
    }

    flist$distinfo <- function(){
        retlist <- bh$distinfo()
        retlist$npars <- retlist$npars - 1
        retlist$parnames <- retlist$parnames[-idx]
        retlist$trans <- log        # assumes this is the trans argument for bh
        retlist$itrans <- exp       # assumes this is the trans argument for bh
        retlist$jacobian <- exp     # assumes this is the trans argument for bh
        retlist$hessian <- retlist$hessian[[-idx]]
        return(retlist)
    }


    flist$basehazard <- function(pars){
        pars <- insert(pars,idx,fixval)
        return(bh$basehazard(pars))
    }

    flist$gradbasehazard <- function(pars){
        pars <- insert(pars,idx,fixval)
        funk <- bh$gradbasehazard(pars) #a function of t
        fun <- function(t,...){
            ans <- funk(t,...) # do the evaluation
            ans <- ans[,-idx]
            return(ans)
        }
        return(fun)
    }

    flist$hessbasehazard <- function(pars){
        pars <- insert(pars,idx,fixval)
        funk <- bh$hessbasehazard(pars) #a function of t
        fun <- function(t,...){
            ans <- funk(t,...) # do the evaluation
            ans <- lapply(ans,function(x){x[-idx,-idx]})
            return(ans)
        }
        return(fun)
    }


    flist$cumbasehazard <- function(pars){
        pars <- insert(pars,idx,fixval)
        return(bh$cumbasehazard(pars))
    }


    flist$gradcumbasehazard <- function(pars){
        pars <- insert(pars,idx,fixval)
        funk <- bh$gradcumbasehazard(pars) #a function of t
        fun <- function(t,...){
            ans <- funk(t,...) # do the evaluation
            ans <- ans[,-idx]
            return(ans)
        }
        return(fun)
    }

    flist$hesscumbasehazard <- function(pars){
        pars <- insert(pars,idx,fixval)
        funk <- bh$hesscumbasehazard(pars) #a function of t
        fun <- function(t,...){
            ans <- funk(t,...) # do the evaluation
            ans <- lapply(ans,function(x){x[-idx,-idx]})
            return(ans)
        }
        return(fun)
    }

    flist$densityquantile <- function(pars,other){
        #fun <- function(probs,...){
            stop("densityquantile not available yet")
            #return((-log(1-probs)/(pars[2]*other$expXbetaplusY))^(1/pars[1]))
        #}
    }


    class(flist) <- c("basehazardspec","list")
    return(flist)

}
