##' BsplineHaz function
##'
##' A function to define a parametric proportional hazards model where the baseline hazard is modelled by a basis spline. 
##' This function returns an object inheriting class 'basehazardspec', list of functions 'distinfo', 'basehazard', 'gradbasehazard', 'hessbasehazard',
##' 'cumbasehazard', 'gradcumbasehazard', 'hesscumbasehazard' and 'densityquantile'
##' 
##' The \code{distinfo} function is used to provide basic distribution specific information to other \code{spatsurv} functions. The user is required 
##' to provide the following information in the returned list: \code{npars}, the number of parameters in this distribution; \code{parnames}, 
##' the names of the parameters; \code{trans}, the transformation scale on which the priors will be provided; \code{itrans}, the inverse 
##' transformation function that will be applied to the parameters before the hazard, and other functions are evaluated; \code{jacobian}, 
##' the derivative of the inverse transformation function with respect to each of the parameters; and \code{hessian}, the second derivatives 
##' of the inverse transformation function with respect to each of the parameters -- note that currently the package \code{spatsurv} 
##' only allows the use of functions where the parameters are transformed independently.
##' 
##' The \code{basehazard} function is used to evaluate the baseline hazard function for the distribution of interest. It returns a 
##' function that accepts as input a vector of times, \code{t} and returns a vector.
##' 
##' The \code{gradbasehazard} function is used to evaluate the gradient of the baseline hazard function with respect to the parameters, 
##' this typically returns a vector. It returns a function that accepts as input a vector of times, \code{t}, and returns a matrix.
##' 
##' The \code{hessbasehazard} function is used to evaluate the Hessian of the baseline hazard function. It returns a function that accepts 
##' as input a vector of times, \code{t} and returns a list of hessian matrices corresponding to each \code{t}.
##' 
##' The \code{cumbasehazard} function is used to evaluate the cumulative baseline hazard function for the distribution of interest. 
##' It returns a function that accepts as input a vector of times, \code{t} and returns a vector.
##' 
##' The \code{gradcumbasehazard} function is used to evaluate the gradient of the cumulative baseline hazard function with respect 
##' to the parameters, this typically returns a vector. It returns a function that accepts as input a vector of times, \code{t}, and returns a matrix.
##' 
##' The \code{hesscumbasehazard} function is used to evaluate the Hessian of the cumulative baseline hazard function. It returns a 
##' function that accepts as input a vector of times, \code{t} and returns a list of hessian matrices corresponding to each \code{t}.
##' 
##' The \code{densityquantile} function is used to return quantiles of the density function. This is NOT REQUIRED for running the MCMC, 
##' merely for us in post-processing with the \code{predict} function where \code{type} is 'densityquantile'. In the case of the Weibull 
##' model for the baseline hazard, it can be shown that the q-th quantile is: 
##'
##' @param times vector of survival times (both censored and uncensored)
##' @param knots vector of knots in ascending order, must include minimum and maximum values of 'times'
##' @param degree degree of the spline basis, default is 3
##' @param MLinits optional starting values for the non-spatial maximisation routine using optim. Note that we are working with the log of the parameters. Default is -10 for each parameter.
##' @return an object inheriting class 'basehazardspec'
##' @seealso \link{exponentialHaz}, \link{gompertzHaz}, \link{makehamHaz}, \link{weibullHaz} 
##' @export

BsplineHaz <- function(times,knots=quantile(times),degree=3,MLinits=NULL){

	polymult <- function(poly1,poly2){
		l1 <- length(poly1)
		l2 <- length(poly2)
		tms <- outer(poly1,poly2)
		ord <- outer(0:(l1-1),0:(l2-1),"+")
		mxord <- ord[l1,l2]
		poly <- c()
		poly <- sapply(0:mxord,function(i){poly[i] <<- sum(tms[ord==i])})
		mxord <- mxord + 1
		if(poly[mxord]==0){
			while(poly[mxord]==0 & length(poly)>1){
				poly <- poly[-mxord]
				mxord <- mxord - 1
			}
		}
		return(poly)
	}

	polyadd <- function(poly1,poly2){
		ans <- rep(0,max(length(poly1),length(poly2)))
		ans[1:length(poly1)] <- ans[1:length(poly1)] + poly1
		ans[1:length(poly2)] <- ans[1:length(poly2)] + poly2
		return(ans)
	}

	alpha <- function(i,j,knots,knotidx){
		if(knots[knotidx==(i+j)]==knots[knotidx==i]){
			return(0)
		}
		else{
			return(c(-knots[knotidx==i]/(knots[knotidx==(i+j)]-knots[knotidx==i]),1/(knots[knotidx==(i+j)]-knots[knotidx==i])))
		}
	}

	B <- function(x,i,j,knots){
		knotidx <- 0:(length(knots)-1)
		if(j==0){
			if(x>= knots[knotidx==i] & x<=knots[knotidx==(i+1)]){
				return(1)
			}
			else{
				return(0)
			}
		}
		else{
			p1 <- alpha(i,j,knots,knotidx)
			p2 <- -alpha(i+1,j,knots,knotidx)
			p2[1] <- p2[1]+1	
			return(polyadd(polymult(p1,B(x,i,j-1,knots)),polymult(p2,B(x,i+1,j-1,knots))))
		}
	}

	midpts <- function(x){
		diffs <- diff(x)
		return(x[1:(length(x)-1)]+diffs/2)
	}

	getBbasis <- function(x,knots,degree){

		if(!all(sort(knots)==knots)){
			stop("Knots must be in ascending order")
		}
		if(knots[1]!=min(x) | rev(knots)[1]!=max(x)){
			stop("First and last knots must be respectively the minimum and maximum of x")
		}

		evalpts <- midpts(knots)
		augknots <- c(rep(knots[1],degree),knots,rep(rev(knots)[1],degree))
		
		polylist <- list()
		for(i in 0:(length(knots)-2+degree)){
			polys <- c()
			for(j in 1:length(evalpts)){
				b <- B(evalpts[j],i,degree,knots=augknots)
				l <- length(b)
				if(l<(degree+1)){
					b <- c(b,rep(0,degree+1-l))
				}
				polys <- rbind(polys,b)			
			}
			rownames(polys) <- NULL
			polylist[[i+1]] <- as.matrix(polys)
		}
		return(list(knots=knots,poly=polylist))
	}

	Bspline.construct <- function(x,basis){
		idx <- as.numeric(cut(x,basis$knots,include.lowest=TRUE))
		xpows <- outer(x,0:(ncol(basis$poly[[1]])-1),"^")
		ans <- t(sapply(1:length(x),function(i){sapply(basis$poly,function(co){sum(co[idx[i],]*xpows[i,])})}))
		return(ans)
	}

	cumulativeBspline.construct <- function(x,basis){

		knots <- basis$knots

		coeffs <- basis$poly

		nb <- length(coeffs)
		np <- ncol(coeffs[[1]])

		addfun <- function(pars){
			cof <- Reduce("+",mapply("*",basis$poly,pars,SIMPLIFY=FALSE))
			if(knots[1]==0){
				return(0)
			}
			else{
				return(knots[1]*sum(cof[1,]*knots[1]^(0:(np-1))))
			}	
		}

		coeffs <- lapply(coeffs,function(mat){return(t(apply(mat,1,function(x){x/(1:np)})))})
		coeffs <- lapply(coeffs,function(mat){return(cbind(0,mat))})
		powers <- 0:np

		idx <- as.numeric(cut(x,knots,include.lowest=TRUE))
		
		idxmax <- max(idx,na.rm=TRUE)

		ints <- t(sapply(coeffs,function(mat){c(0,sapply(2:idxmax,function(i){sum(mat[i-1,]*(knots[i]^powers-knots[i-1]^powers))}))}))
		cumints <- t(apply(ints,1,cumsum))

		integ <- c()
		for (i in 1:length(idx)){
			xpow <- x[i]^powers
			kpow <- knots[idx[i]]^powers
			integ <- rbind(integ,cumints[,idx[i]] + sapply(1:length(coeffs),function(j){sum(coeffs[[j]][idx[i],]*(xpow-kpow))}))
		}
		
		return(list(integral=integ,toadd=addfun))
	}

	basis <- getBbasis(x=times,knots,degree=degree)

	tpows <- outer(times,0:degree,"^")

	basismatrix <- Bspline.construct(x=times,basis=basis) 
	basismatrix[basismatrix<0] <- 0 # these are small negative numbers

	cbs <- cumulativeBspline.construct(x=times,basis=basis) 
	
	#if(any(times<basis$knots[1])){
	#	stop("All 'times' must lie between smallest and largest knot values")
	#}

	np <- length(basis$poly)

    flist <- list()

    flist$distinfo <- function(){
        retlist <- list()
        retlist$npars <- np
        retlist$parnames <- paste("lambda",1:np,sep="")
        retlist$trans <- log
        retlist$itrans <- exp
        retlist$jacobian <- exp
        retlist$hessian <- lapply(1:np,function(zz){return(exp)})
        if(is.null(MLinits)){
        	retlist$MLinits <- rep(-10,np)	
        }
        else{
        	retlist$MLinits <- MLinits
        }        
        return(retlist)
    }

    test <- flist$distinfo()
    cat("Using B-spline with ",test$npars," parameters.\n")
    
    flist$basehazard <- function(pars){
        fun <- function(t){
        	idx <- match(t,times)
        	if(any(is.na(idx))){
        		basismatrix <- Bspline.construct(x=t,basis=basis) 
        		basismatrix[basismatrix<0] <- 0 # these are small negative numbers
        		return(colSums(pars*t(basismatrix)))
        	}
        	else{
        		return(colSums(pars*t(basismatrix[idx,])))
        	}
        }
        return(fun)  
    }
    
    flist$gradbasehazard <- function(pars){
        fun <- function(t){
        	idx <- match(t,times)
        	return(basismatrix[idx,])
        }
        return(fun)        
    }
    
    flist$hessbasehazard <- function(pars){
        funfun <- function(t,pars){ 
            return(matrix(0,np,np)) 
        }
        
        fun <- function(t){
            return(lapply(t,funfun,pars=pars))
        }
        return(fun)
        
    }
    
    flist$cumbasehazard <- function(pars){
        fun <- function(t){
        	idx <- match(t,times)       	
        	if(any(is.na(idx))){
        		cbs <- cumulativeBspline.construct(x=t,basis=basis) 
        		return(colSums(pars*t(cbs$integral)) + cbs$toadd(pars))
        	}
        	else{
        		return(colSums(pars*t(cbs$integral[idx,])) + cbs$toadd(pars))
        	}
        }
        return(fun) 
    }
    
    flist$gradcumbasehazard <- function(pars){
        fun <- function(t){
            idx <- match(t,times)
        	return(cbs$integral[idx,])        
        }
        return(fun)    
    }
    
    flist$hesscumbasehazard <- function(pars){
        funfun <- function(t,pars){
            return(matrix(0,np,np)) 
        }
        
        fun <- function(t){
            return(lapply(t,funfun,pars=pars))
        }
        return(fun)
    }
    
    flist$densityquantile <- function(pars,other){
        fun <- function(probs){
            stop("densityquantile not available yet")
            #return((-log(1-probs)/(pars[2]*other$expXbetaplusY))^(1/pars[1]))
        }
        return(fun)    
    }


    class(flist) <- c("basehazardspec","list")
    return(flist)
}