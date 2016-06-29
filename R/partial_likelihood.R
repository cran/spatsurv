##' timevaryingPL function
##'
##' A function to 
##'
##' @param formula a formula of the form 'S ~ coef1 + coef2' etc the object S will be created 
##' @param t0 X
##' @param t X
##' @param delta censoring indicator a vector of 1 for an event and 0 for censoring
##' @param dist X 
##' @param data X 
##' @param ties X default is Efron
##' @param optimcontrol X 
##' @return ...
##' @export

timevaryingPL <- function(formula,t0,t,delta,dist,data,ties="Efron",optimcontrol=NULL){

    interv <- interval(t0,t0+t)

    # n <- 100
    # rg <- range(as.numeric(c(t0,t0+t)))
    # plot(NULL,xlim=rg,ylim=c(0,n+1))
    # for(i in 1:n){
    #     lines(c(as.numeric(t0[i]),as.numeric(t0[i]+t[i])),c(i,i))
    # }

    


    data$S <- 1 # dummy

    # ##########
    # # This chunk of code borrowed from flexsurvreg    
    # ##########  
    # call <- match.call()
    # indx <- match(c("formula", "data"), names(call), nomatch = 0)
    # if (indx[1] == 0){ 
    #     stop("A \"formula\" argument is required")
    # }
    # temp <- call[c(1, indx)]
    # temp[[1]] <- as.name("model.frame")
    # m <- eval(temp, parent.frame())

    # Terms <- attr(m, "terms")
    # Z <- model.matrix(Terms, m)
    # ##########
    # # End of borrowed code    
    # ##########

    Z <- model.matrix(formula,data)
    Z <- Z[, -1, drop = FALSE] 

    if(nrow(Z)!=length(interv)){
        stop("missing covariate data not allowed")
    }

    nb <- ncol(Z)
    no <- distinfo(dist)()$npars
    np <- nb + no

    n <- nrow(Z)

    Period2numeric <- function(x){
        return(day(x) + hour(x)/24 + minute(x)/(24*60) + second(x)/(24*60*60))
    }

    plfun <- function(pars){
        print(pars)
        beta <- pars[1:nb]
        omega <- pars[(nb+1):np]
        #omegaorig <- omega # recall we are working with omega on the transformed scale
        omega <- distinfo(dist)()$itrans(omega) # this is omega on the correct scale
        haz <- setupHazard(dist=dist,pars=omega,grad=FALSE,hess=FALSE)
        Zbeta <- Z%*%beta
        expZbeta <- exp(Zbeta)

        t0_plus_t <- t0+t
        dec_t <- Period2numeric(t)

        h0t <- haz$h(dec_t)
        logh0t <- log(h0t)

        logPL <- 0
        for(i in 1:n){
            if(delta[i]==0){
                next
            }
            else{
                contrib <- logh0t[i] + Zbeta[i]
                RS <- t0_plus_t[i]%within%interv
                #print(sum(RS))
                tij <- as.numeric(difftime(t0_plus_t[i],t0[RS]))
                contrib <- contrib - log(sum(haz$h(tij)*expZbeta[RS]))
                logPL <- logPL + contrib
            }
        }

        #browser()

        #print(logPL)

        return(logPL)

    }

    plfundrv <- function(pars){
        beta <- pars[1:nb]
        omega <- pars[(nb+1):np]
        #omegaorig <- omega # recall we are working with omega on the transformed scale
        omega <- distinfo(dist)()$itrans(omega) # this is omega on the correct scale
        haz <- setupHazard(dist=dist,pars=omega,grad=TRUE,hess=FALSE)
        Zbeta <- Z%*%beta
        expZbeta <- exp(Zbeta)

        t0_plus_t <- t0 + t
        dec_t <- Period2numeric(t)

        h0t <- haz$h(dec_t)
        gradh0t <- haz$gradh(dec_t)
        hquo <- gradh0t / h0t
        #logh0t <- log(h0t)

        dl_dbeta <- rep(0,nb)
        dl_domega <- rep(0,no)
        for(i in 1:n){
            if(delta[i]==0){
                next
            }
            else{
                RS <- t0_plus_t[i]%within%interv
                tij <- as.numeric(difftime(t0_plus_t[i],t0[RS]))
                htij <- haz$h(tij)
                gradhtij <- haz$gradh(tij)
                fij <- htij*expZbeta[RS]
                dfij_dbeta <- fij*Z[RS,] #htij*expZbeta[RS]*Z[RS,]
                dfij_domega <- expZbeta[RS]*gradhtij
                
                if(nb>1){
                    dl_dbeta <- dl_dbeta + Z[i,] - colSums(dfij_dbeta) / sum(fij)
                }
                else{
                    dl_dbeta <- dl_dbeta + Z[i] - sum(dfij_dbeta) / sum(fij)
                }

                if(no>1){
                    dl_domega <- dl_domega + hquo[i,] - colSums(dfij_domega) / sum(fij)
                }
                else{
                    dl_domega <- dl_domega + hquo[i] - sum(dfij_domega) / sum(fij)
                }

                #browser()
            }
        }

        #browser()

        return(c(dl_dbeta,dl_domega)) 

    }

    opt <- optim(rep(0,np),fn=plfun,gr=plfundrv,method="BFGS",control=list(fnscale=-1))
    #opt1 <- optim(rep(0,np),fn=plfun,control=list(fnscale=-1))
    

    browser()

    
}