##' TwoWayHazAdditive function
##'
##' A function to
##'
##' @param bhlist X
##' @param bhtime X
##' @param bhfix X
##' @param MLinits X
##' @return ...
##' @export

TwoWayHazAdditive <- function(bhlist,bhtime,bhfix,MLinits=NULL){

    stop("I think there is an issue with the TwoWayHazAdditive function ...")

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
            ans <- sapply(1:nhaz,function(i){bhlist[[i]]$basehazard(pars[paridx[[i]]])(bhtime[[i]])})
            return(apply(ans,1,sum))
        }
        return(fun)
    }

    flist$gradbasehazard <- function(pars){
        fun <- function(t,...){
            hdash <- sapply(1:nhaz,function(i){bhlist[[i]]$gradbasehazard(pars[paridx[[i]]])(bhtime[[i]])})
            return(cbind(hdash[[1]],hdash[[2]]))
        }
        return(fun)
    }

    flist$hessbasehazard <- function(pars){
        fun <- function(t,...){
            hess <- lapply(1:nhaz,function(i){bhlist[[i]]$hessbasehazard(pars[paridx[[i]]])(bhtime[[i]])}) # here and above 2 lines, i is bh index
            getMlist <- function(i){ # here i is subject index
                return(lapply(1:nhaz,function(j){hess[[j]][[i]]}))
            }
            hesslist <- list()
            for(i in 1:length(bhtime[[1]])){ # loop over subjects, yes this is likely slow, but does not matter ...
                m <- blockDiag(getMlist(i))
                hesslist[[i]] <- m
            }
            return(hesslist)
        }
        return(fun)

    }

    flist$cumbasehazard <- function(pars){
        fun <- function(t,...){
            tms <- sapply(bhtime,cbind)
            H <- sapply(1:nhaz,function(i){bhlist[[i]]$cumbasehazard(pars[paridx[[i]]])(bhtime[[i]])})
            return(tms[,1]*H[,2] + tms[,2]*H[,1])
        }
        return(fun)
    }

    flist$gradcumbasehazard <- function(pars){
        fun <- function(t,...){
            tms <- sapply(bhtime,cbind)
            H <- sapply(1:nhaz,function(i){bhlist[[i]]$cumbasehazard(pars[paridx[[i]]])(bhtime[[i]])})
            Hdash <- sapply(1:nhaz,function(i){bhlist[[i]]$gradcumbasehazard(pars[paridx[[i]]])(bhtime[[i]])})
            return(cbind(H[,2]+tms[,2]*Hdash[[1]] , H[,1]+tms[,1]*Hdash[[2]]))
        }
        return(fun)
    }

    flist$hesscumbasehazard <- function(pars){
        fun <- function(t,...){
            hess <- lapply(1:nhaz,function(i){bhlist[[i]]$hesscumbasehazard(pars[paridx[[i]]])(bhtime[[i]])}) # here and above 2 lines, i is bh index
            hess[[1]] <- mapply("*",hess[[1]],bhtime[[2]],SIMPLIFY=FALSE)
            hess[[2]] <- mapply("*",hess[[2]],bhtime[[1]],SIMPLIFY=FALSE)
            Hdash <- sapply(1:nhaz,function(i){as.matrix(bhlist[[i]]$gradcumbasehazard(pars[paridx[[i]]])(bhtime[[i]]))})
            getMlist <- function(i){ # here i is subject index
                return(lapply(1:nhaz,function(j){hess[[j]][[i]]}))
            }
            hesslist <- list()
            for(i in 1:length(bhtime[[1]])){ # loop over subjects, yes this is likely slow, but does not matter ...
                m <- blockDiag(getMlist(i))
                od <- outer(Hdash[[1]][i,],Hdash[[2]][i,],"+")
                m[3,1:2] <- m[1:2,3] <- od
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
