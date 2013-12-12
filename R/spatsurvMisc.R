##' gensens function
##'
##' A function to generate observed times given a vector of true survival times and a vector of censoring times. Used in the simulation of
##' survival data
##'
##' @param survtimes a vector of survival times 
##' @param censtimes a vector of censoring times 
##' @return a named list containing 'obstimes', the observed time of the event; and 'censored', the censoring indicator which is equal to 1 if the
##' event is observed and 0 otherwise.
##' @export


gensens <- function(survtimes,censtimes){
    
    n <- length(survtimes)
    
    if(length(survtimes)!=length(censtimes)){
        stop("survtimes and censtimes should have the same length")
    } 

    obstimes <- survtimes
    cens <- rep(1,n)
   
    for(i in 1:n){
        if(censtimes[i]<survtimes[i]){
            obstimes[i] <- censtimes[i]
            cens[i] <- 0
        }
    }
    
    return(Surv(time=obstimes,event=cens))
}






##' plotsurv function
##'
##' A function to produce a 2-D plot of right censored spatial survival data.
##'
##' @param spp A spatial points data frame
##' @param ss A Surv object (with right-censoring) 
##' @param maxcex maximum size of dots default is equavalent to setting cex equal to 1
##' @param background a background object to plot default is null, which gives a blamk background note that if non-null, the parameters xlim and ylim will be derived from this object.
##' @param eventpt The type of point to illustrate events, default is 19 (see ?pch) 
##' @param eventcol the colour of events, default is black
##' @param censpt The type of point to illustrate events, default is "+" (see ?pch)
##' @param censcol the colour of censored observations, default is red
##' @param xlim optional x-limits of plot, default is to choose this automatically 
##' @param ylim optional y-limits of plot, default is to choose this automatically

##' @param ... other arguments to pass to plot
##' @return Plots the survival data non-censored observations appear as dots and censored observations as crosses. The size of the dot is proportional to the observed time.
##' @export

plotsurv <- function(spp,ss,maxcex=1,background=NULL,eventpt=19,eventcol="red",censpt="+",censcol="black",xlim=NULL,ylim=NULL,...){
    crds <- coordinates(spp)
    if(is.null(xlim)){
        if(is.null(background)){
            xlim <- range(crds[,1])
        }
    }
    if(is.null(ylim)){
        if(is.null(background)){
            ylim <- range(crds[,2])
        }
    }

    event <- ss[,"status"] == 1 # event indicator
    cexx <- maxcex* ss[,"time"] / max(ss[,"time"])    
    
    plot(background,xlim=xlim,ylim=ylim,...)
    points(crds[event,],pch=eventpt,col=eventcol,cex=cexx[event])
    points(crds[!event,],pch=censpt,col=censcol,cex=cexx[!event])
    
}


##' inference.control function
##'
##' A function to control inferential settings. 
##'
##' @param gridded logical. Whether to perform compuation on a grid. Default is FALSE. Note in version 0.9-1, this is still in the testing phase.
##' @param cellwidth the width of computational cells to use 
##' @param ext integer the number of times to extend the computational grid by in order to perform compuitation. The default is 2. 
##' @return returns parameters to be used in the function survspat
##' @seealso \link{survspat}
##' @export

inference.control <- function(gridded=FALSE,cellwidth=NULL,ext=2){
    ans <- list()
    ans$gridded <- gridded
    ans$cellwidth <- cellwidth 
    ans$ext <- ext 
    class(ans) <- c("inference.control","list")
    return(ans)
}