##' baselinehazard_multiWay function
##'
##' A function to
##'
##' @param x X
##' @param probs X
##' @param cumulative X
##' @param plot X
##' @param joint X
##' @param xlims X
##' @param ylims X
##' @param ... X
##' @return ...
##' @export


baselinehazard_multiWay <- function(x,probs=c(0.025,0.5,0.975),cumulative=FALSE,plot=TRUE,joint=FALSE,xlims=NULL,ylims=NULL,...){

    if(length(probs)!=3){
        stop("length(probs)!=3 ... oh no you don't :-)")
    }

    omegasamp <- x$omegasamp
    env <- environment(x$dist$distinfo)
    assign("post_processing_mode",TRUE,envir=env)
    nhaz <- get("nhaz",envir=env)
    times <- get("bhtime",envir=env)
    tmsord <- lapply(times,order)
    nobs <- length(times[[1]])

    if(nhaz!=2){
        stop("nhaz!=2 ... no can do :-)")
    }

    fun <- function(pars){
        f <- basehazard(x$dist)(pars)
        return(f(x$survivaldata[,"time"]))
    }
    YLAB <- "Baseline Hazard"

    if(cumulative){
        fun <- function(pars){
            f <- cumbasehazard(x$dist)(pars)
            return(f(x$survivaldata[,"time"]))
        }
        YLAB <- "Cumulative Baseline Hazard"
    }

    samp <- t(apply(omegasamp,1,fun))

    #browser()
    toreturn <- list()
    for(i in 1:nhaz){
        toreturn[[i]] <- t(apply(samp[,((i-1)*nobs+1):(i*nobs)],2,quantile,probs=probs,na.rm=TRUE))
    }
    toreturn[[nhaz+1]] <- times

    # now reorder in order to get to plot
    for(i in 1:nhaz){
        toreturn[[i]] <- toreturn[[i]][tmsord[[i]],]
        toreturn[[nhaz+1]][[i]] <- toreturn[[nhaz+1]][[i]][tmsord[[i]]]
    }


    names(toreturn) <- c(paste("hazard",1:nhaz,sep=""),"times")
    if(plot){
        if(joint){
            # laymat <- matrix(1,nhaz,nhaz+2)
            # laymat[1,(nhaz+1):(nhaz+2)] <- 2
            # laymat[2,(nhaz+1):(nhaz+2)] <- 3
            # layout(laymat)
            image(toreturn[[nhaz+1]][[1]],toreturn[[nhaz+1]][[2]],outer(toreturn[[1]][,2],toreturn[[2]][,2],FUN="*"),xlab="time scale 1",ylab="time scale 2")
            for(i in 1:nhaz){
                matplot(toreturn[[nhaz+1]][[i]],toreturn[[i]],type="l",col=c("black","black","black"),lty=c("dotted","solid","dashed"),xlab=paste("time scale",i),ylab=YLAB,xlim=xlims[[i]],ylim=ylims[[i]],...)
                legend("topright",lty=c("dashed","solid","dotted"),col=rev(c("black","black","black")),legend=rev(probs))
            }
        }
        else{
            #layout(matrix(c(1,1,1,1,2,2,2,2),2,4))
            #par(mfrow=c(1,2))
            for(i in 1:nhaz){
                matplot(toreturn[[nhaz+1]][[i]],toreturn[[i]],type="l",col=c("black","black","black"),lty=c("dotted","solid","dashed"),xlab=paste("time scale",i),ylab=YLAB,xlim=xlims[[i]],ylim=ylims[[i]],...)
                legend("topright",lty=c("dashed","solid","dotted"),col=rev(c("black","black","black")),legend=rev(probs))
            }
        }
    }


    return(list(t=t,qts=toreturn))
}
