# ##' TPSHaz function
# ##'
# ##' A function to define a parametric proportional hazards model where the two-way baseline hazard is modelled by a tensor product spline.
# ##' This function returns an object inheriting class 'basehazardspec', list of functions 'distinfo', 'basehazard', 'gradbasehazard', 'hessbasehazard',
# ##' 'cumbasehazard', 'gradcumbasehazard', 'hesscumbasehazard' and 'densityquantile'
# ##'
# ##' The \code{distinfo} function is used to provide basic distribution specific information to other \code{spatsurv} functions. The user is required
# ##' to provide the following information in the returned list: \code{npars}, the number of parameters in this distribution; \code{parnames},
# ##' the names of the parameters; \code{trans}, the transformation scale on which the priors will be provided; \code{itrans}, the inverse
# ##' transformation function that will be applied to the parameters before the hazard, and other functions are evaluated; \code{jacobian},
# ##' the derivative of the inverse transformation function with respect to each of the parameters; and \code{hessian}, the second derivatives
# ##' of the inverse transformation function with respect to each of the parameters -- note that currently the package \code{spatsurv}
# ##' only allows the use of functions where the parameters are transformed independently.
# ##'
# ##' The \code{basehazard} function is used to evaluate the baseline hazard function for the distribution of interest. It returns a
# ##' function that accepts as input a vector of times, \code{t} and returns a vector.
# ##'
# ##' The \code{gradbasehazard} function is used to evaluate the gradient of the baseline hazard function with respect to the parameters,
# ##' this typically returns a vector. It returns a function that accepts as input a vector of times, \code{t}, and returns a matrix.
# ##'
# ##' The \code{hessbasehazard} function is used to evaluate the Hessian of the baseline hazard function. It returns a function that accepts
# ##' as input a vector of times, \code{t} and returns a list of hessian matrices corresponding to each \code{t}.
# ##'
# ##' The \code{cumbasehazard} function is used to evaluate the cumulative baseline hazard function for the distribution of interest.
# ##' It returns a function that accepts as input a vector of times, \code{t} and returns a vector.
# ##'
# ##' The \code{gradcumbasehazard} function is used to evaluate the gradient of the cumulative baseline hazard function with respect
# ##' to the parameters, this typically returns a vector. It returns a function that accepts as input a vector of times, \code{t}, and returns a matrix.
# ##'
# ##' The \code{hesscumbasehazard} function is used to evaluate the Hessian of the cumulative baseline hazard function. It returns a
# ##' function that accepts as input a vector of times, \code{t} and returns a list of hessian matrices corresponding to each \code{t}.
# ##'
# ##' The \code{densityquantile} function is used to return quantiles of the density function. This is NOT REQUIRED for running the MCMC,
# ##' merely for us in post-processing with the \code{predict} function where \code{type} is 'densityquantile'. In the case of the Weibull
# ##' model for the baseline hazard, it can be shown that the q-th quantile is:
# ##'
# ##' @param times1 vector of survival times (both censored and uncensored)
# ##' @param knots vector of knots in ascending order, must include minimum and maximum values of 'times'
# ##' @param degree degree of the spline basis, default is 3
# ##' @param MLinits optional starting values for the non-spatial maximisation routine using optim. Note that we are working with the log of the parameters. Default is -10 for each parameter.
# ##' @return an object inheriting class 'basehazardspec'
# ##' @seealso \link{exponentialHaz}, \link{gompertzHaz}, \link{makehamHaz}, \link{weibullHaz}
# ##' @export
#
# TPSHaz <- function(times1,times2,knots1=quantile(times1),knots2=quantile(times2),degree1=3,degree2=3,MLinits=NULL){
#
#     knots1[1] <- 0 # helps with definition of cumulative splines
#     knots2[1] <- 0
#
# 	basis1 <- getBbasis(x=times1,knots1,degree=degree1,force=TRUE)
#     basis2 <- getBbasis(x=times2,knots2,degree=degree2,force=TRUE)
#
# 	basismatrix1 <- Bspline.construct(x=times1,basis=basis1)
# 	basismatrix1[basismatrix1<0] <- 0 # these are small negative numbers
#     basismatrix2 <- Bspline.construct(x=times2,basis=basis2)
# 	basismatrix2[basismatrix2<0] <- 0 # these are small negative numbers
#
# 	cbs1 <- cumulativeBspline.construct(x=times1,basis=basis1)
#     cbs2 <- cumulativeBspline.construct(x=times2,basis=basis2)
#
# 	np1 <- length(basis1$poly)
#     np2 <- length(basis2$poly)
#
#     flist <- list()
#
#     flist$distinfo <- function(){
#         retlist <- list()
#         retlist$npars <- np1
#         retlist$parnames <- paste("lambda",1:(np1),sep="")
#         retlist$trans <- log
#         retlist$itrans <- exp
#         retlist$jacobian <- exp
#         retlist$hessian <- lapply(1:(np1),function(zz){return(exp)})
#         if(is.null(MLinits)){
#         	retlist$MLinits <- rep(-10,np1)
#         }
#         else{
#         	retlist$MLinits <- MLinits
#         }
#         return(retlist)
#     }
#
#     test <- flist$distinfo()
#     cat(paste("Using Tensor Product Spline with ",np1,"x",np2,"=",test$npars," parameters.\n",sep=""))
#
#     flist$basehazard <- function(pars){
#         fun <- function(t,...){
#         	#idx <- match(t,times) # NOT NECESSARY (?)
#         	#if(any(is.na(idx))){
#         	#	basismatrix <- Bspline.construct(x=t,basis=basis)
#         	#	basismatrix[basismatrix<0] <- 0 # these are small negative numbers
#         	#	return(colSums(pars*t(basismatrix)))
#         	#}
#         	#else{
#                 M <- diag(pars)#matrix(pars,np1,np2)
#                 #browser()
#         		return(apply(cbind(basismatrix1,basismatrix2),1,function(x){x[1:np1]%*%(M%*%x[(np1+1):(np1+np2)])}))
#         	#}
#         }
#         return(fun)
#     }
#
#     flist$gradbasehazard <- function(pars){
#         fun <- function(t,...){
#         	#idx <- match(t,times)
#         	return(apply(cbind(basismatrix1,basismatrix2),1,function(x){as.vector(diag(outer(x[1:np1],x[(np1+1):(np1+np2)])))}))
#         }
#         return(fun)
#     }
#
#     flist$hessbasehazard <- function(pars){
#         funfun <- function(t,pars){
#             return(matrix(0,np1,np1))
#         }
#
#         fun <- function(t,...){
#             return(lapply(t,funfun,pars=pars))
#         }
#         return(fun)
#
#     }
#
#     flist$cumbasehazard <- function(pars){
#         fun <- function(t,...){
#         	#idx <- match(t,times)
#         	#if(any(is.na(idx))){
#         	#	cbs <- cumulativeBspline.construct(x=t,basis=basis)
#         	#	return(colSums(pars*t(cbs$integral)) + cbs$toadd(pars))
#         	#}
#         	#else{
#                 M <- diag(pars)#matrix(pars,np1,np2)
#                 #browser()
#         		return(apply(cbind(cbs1$integral,cbs2$integral),1,function(x){x[1:np1]%*%(M%*%x[(np1+1):(np1+np2)])})) #  + cbs1$toadd(pars)*cbs2$toadd(pars) NOT REQUIRED AS KNOTS[1] = 0        	#}
#         }
#         return(fun)
#     }
#
#     flist$gradcumbasehazard <- function(pars){
#         fun <- function(t,...){
#             #idx <- match(t,times)
#         	#return(cbs$integral[idx,])
#             return(apply(cbind(cbs1$integral,cbs1$integral),1,function(x){as.vector(diag(outer(x[1:np1],x[(np1+1):(np1+np2)])))}))
#         }
#         return(fun)
#     }
#
#     flist$hesscumbasehazard <- function(pars){
#         funfun <- function(t,pars){
#             return(matrix(0,np1,np1))
#         }
#
#         fun <- function(t,...){
#             return(lapply(t,funfun,pars=pars))
#         }
#         return(fun)
#     }
#
#     flist$densityquantile <- function(pars,other){
#         fun <- function(probs,...){
#             stop("densityquantile not available yet")
#             #return((-log(1-probs)/(pars[2]*other$expXbetaplusY))^(1/pars[1]))
#         }
#         return(fun)
#     }
#
#
#     class(flist) <- c("basehazardspec","list")
#     return(flist)
# }
#
#
# ##' TPSHaz1 function
# ##'
# ##' A function to define a parametric proportional hazards model where the two-way baseline hazard is modelled by a tensor product spline.
# ##' This function returns an object inheriting class 'basehazardspec', list of functions 'distinfo', 'basehazard', 'gradbasehazard', 'hessbasehazard',
# ##' 'cumbasehazard', 'gradcumbasehazard', 'hesscumbasehazard' and 'densityquantile'
# ##'
# ##' The \code{distinfo} function is used to provide basic distribution specific information to other \code{spatsurv} functions. The user is required
# ##' to provide the following information in the returned list: \code{npars}, the number of parameters in this distribution; \code{parnames},
# ##' the names of the parameters; \code{trans}, the transformation scale on which the priors will be provided; \code{itrans}, the inverse
# ##' transformation function that will be applied to the parameters before the hazard, and other functions are evaluated; \code{jacobian},
# ##' the derivative of the inverse transformation function with respect to each of the parameters; and \code{hessian}, the second derivatives
# ##' of the inverse transformation function with respect to each of the parameters -- note that currently the package \code{spatsurv}
# ##' only allows the use of functions where the parameters are transformed independently.
# ##'
# ##' The \code{basehazard} function is used to evaluate the baseline hazard function for the distribution of interest. It returns a
# ##' function that accepts as input a vector of times, \code{t} and returns a vector.
# ##'
# ##' The \code{gradbasehazard} function is used to evaluate the gradient of the baseline hazard function with respect to the parameters,
# ##' this typically returns a vector. It returns a function that accepts as input a vector of times, \code{t}, and returns a matrix.
# ##'
# ##' The \code{hessbasehazard} function is used to evaluate the Hessian of the baseline hazard function. It returns a function that accepts
# ##' as input a vector of times, \code{t} and returns a list of hessian matrices corresponding to each \code{t}.
# ##'
# ##' The \code{cumbasehazard} function is used to evaluate the cumulative baseline hazard function for the distribution of interest.
# ##' It returns a function that accepts as input a vector of times, \code{t} and returns a vector.
# ##'
# ##' The \code{gradcumbasehazard} function is used to evaluate the gradient of the cumulative baseline hazard function with respect
# ##' to the parameters, this typically returns a vector. It returns a function that accepts as input a vector of times, \code{t}, and returns a matrix.
# ##'
# ##' The \code{hesscumbasehazard} function is used to evaluate the Hessian of the cumulative baseline hazard function. It returns a
# ##' function that accepts as input a vector of times, \code{t} and returns a list of hessian matrices corresponding to each \code{t}.
# ##'
# ##' The \code{densityquantile} function is used to return quantiles of the density function. This is NOT REQUIRED for running the MCMC,
# ##' merely for us in post-processing with the \code{predict} function where \code{type} is 'densityquantile'. In the case of the Weibull
# ##' model for the baseline hazard, it can be shown that the q-th quantile is:
# ##'
# ##' @param times vector of survival times (both censored and uncensored)
# ##' @param knots vector of knots in ascending order, must include minimum and maximum values of 'times'
# ##' @param degree degree of the spline basis, default is 3
# ##' @param MLinits optional starting values for the non-spatial maximisation routine using optim. Note that we are working with the log of the parameters. Default is -10 for each parameter.
# ##' @return an object inheriting class 'basehazardspec'
# ##' @seealso \link{exponentialHaz}, \link{gompertzHaz}, \link{makehamHaz}, \link{weibullHaz}
# ##' @export
#
# TPSHaz1 <- function(times1,times2,knots1=quantile(times1),knots2=quantile(times2),degree1=3,degree2=3,MLinits=NULL){
#
#     knots1[1] <- 0 # helps with definition of cumulative splines
#     knots2[1] <- 0
#
# 	basis1 <- getBbasis(x=times1,knots1,degree=degree1,force=TRUE)
#     basis2 <- getBbasis(x=times2,knots2,degree=degree2,force=TRUE)
#
# 	basismatrix1 <- Bspline.construct(x=times1,basis=basis1)
# 	basismatrix1[basismatrix1<0] <- 0 # these are small negative numbers
#     basismatrix2 <- Bspline.construct(x=times2,basis=basis2)
# 	basismatrix2[basismatrix2<0] <- 0 # these are small negative numbers
#
# 	cbs1 <- cumulativeBspline.construct(x=times1,basis=basis1)
#     cbs2 <- cumulativeBspline.construct(x=times2,basis=basis2)
#
# 	np1 <- length(basis1$poly)
#     np2 <- length(basis2$poly)
#
#     flist <- list()
#
#     flist$distinfo <- function(){
#         retlist <- list()
#         retlist$npars <- np1*np2
#         retlist$parnames <- paste("lambda",1:(np1*np2),sep="")
#         retlist$trans <- log
#         retlist$itrans <- exp
#         retlist$jacobian <- exp
#         retlist$hessian <- lapply(1:(np1*np2),function(zz){return(exp)})
#         if(is.null(MLinits)){
#         	retlist$MLinits <- -runif(np1*np2)# rep(-10,np1*np2)
#         }
#         else{
#         	retlist$MLinits <- MLinits
#         }
#         return(retlist)
#     }
#
#     test <- flist$distinfo()
#     cat(paste("Using Tensor Product Spline with ",np1,"x",np2,"=",test$npars," parameters.\n",sep=""))
#
#     flist$basehazard <- function(pars){
#         fun <- function(t,...){
#         	#idx <- match(t,times) # NOT NECESSARY (?)
#         	#if(any(is.na(idx))){
#         	#	basismatrix <- Bspline.construct(x=t,basis=basis)
#         	#	basismatrix[basismatrix<0] <- 0 # these are small negative numbers
#         	#	return(colSums(pars*t(basismatrix)))
#         	#}
#         	#else{
#                 M <- matrix(pars,np1,np2)
#                 #browser()
#         		return(apply(cbind(basismatrix1,basismatrix2),1,function(x){x[1:np1]%*%(M%*%x[(np1+1):(np1+np2)])}))
#         	#}
#         }
#         return(fun)
#     }
#
#     flist$gradbasehazard <- function(pars){
#         fun <- function(t,...){
#         	#idx <- match(t,times)
#         	return(apply(cbind(basismatrix1,basismatrix2),1,function(x){as.vector(outer(x[1:np1],x[(np1+1):(np1+np2)]))}))
#         }
#         return(fun)
#     }
#
#     flist$hessbasehazard <- function(pars){
#         funfun <- function(t,pars){
#             return(matrix(0,np1*np2,np1*np2))
#         }
#
#         fun <- function(t,...){
#             return(lapply(t,funfun,pars=pars))
#         }
#         return(fun)
#
#     }
#
#     flist$cumbasehazard <- function(pars){
#         fun <- function(t,...){
#         	#idx <- match(t,times)
#         	#if(any(is.na(idx))){
#         	#	cbs <- cumulativeBspline.construct(x=t,basis=basis)
#         	#	return(colSums(pars*t(cbs$integral)) + cbs$toadd(pars))
#         	#}
#         	#else{
#                 M <- matrix(pars,np1,np2)
#                 #browser()
#         		return(apply(cbind(cbs1$integral,cbs2$integral),1,function(x){x[1:np1]%*%(M%*%x[(np1+1):(np1+np2)])})) #  + cbs1$toadd(pars)*cbs2$toadd(pars) NOT REQUIRED AS KNOTS[1] = 0        	#}
#         }
#         return(fun)
#     }
#
#     flist$gradcumbasehazard <- function(pars){
#         fun <- function(t,...){
#             #idx <- match(t,times)
#         	#return(cbs$integral[idx,])
#             return(apply(cbind(cbs1$integral,cbs1$integral),1,function(x){as.vector(outer(x[1:np1],x[(np1+1):(np1+np2)]))}))
#         }
#         return(fun)
#     }
#
#     flist$hesscumbasehazard <- function(pars){
#         funfun <- function(t,pars){
#             return(matrix(0,np1*np2,np1*np2))
#         }
#
#         fun <- function(t,...){
#             return(lapply(t,funfun,pars=pars))
#         }
#         return(fun)
#     }
#
#     flist$densityquantile <- function(pars,other){
#         fun <- function(probs,...){
#             stop("densityquantile not available yet")
#             #return((-log(1-probs)/(pars[2]*other$expXbetaplusY))^(1/pars[1]))
#         }
#         return(fun)
#     }
#
#
#     class(flist) <- c("basehazardspec","list")
#     return(flist)
# }
