##' polymult function
##'
##' A function to multiply two polynomials in the form of vectors of coefficients. The first
##' element of the vector being the constant (order 0) term
##'
##' @param poly1 a vector of coefficients for the first polynomial of length degree plus 1
##' @param poly2 a vector of coefficients for the second polynomial of length degree plus 1
##' @return the coefficients of the product of poly1 and poly2
##' @export

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



##' polyadd function
##'
##' A function to add two polynomials in the form of vectors of coefficients. The first
##' element of the vector being the constant (order 0) term
##'
##' @param poly1 a vector of coefficients for the first polynomial of length degree plus 1
##' @param poly2 a vector of coefficients for the second polynomial of length degree plus 1
##' @return the coefficients of the sum of poly1 and poly2
##' @export

polyadd <- function(poly1,poly2){
    ans <- rep(0,max(length(poly1),length(poly2)))
    ans[1:length(poly1)] <- ans[1:length(poly1)] + poly1
    ans[1:length(poly2)] <- ans[1:length(poly2)] + poly2
    return(ans)
}



##' alpha function
##'
##' A function used in calculating the coefficients of a B-spline curve
##'
##' @param i index i
##' @param j index j
##' @param knots knot vector
##' @param knotidx knot index
##' @return a vector
##' @export

alpha <- function(i,j,knots,knotidx){
    if(knots[knotidx==(i+j)]==knots[knotidx==i]){
        return(0)
    }
    else{
        return(c(-knots[knotidx==i]/(knots[knotidx==(i+j)]-knots[knotidx==i]),1/(knots[knotidx==(i+j)]-knots[knotidx==i])))
    }
}



##' B function
##'
##' A recursive function used in calculating the coefficients of a B-spline curve
##'
##' @param x locations at which to evaluate the B-spline
##' @param i index i
##' @param j index j
##' @param knots a knot vector
##' @return a vector of polynomial coefficients
##' @export

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



##' midpts function
##'
##' A function to compute the midpoints of a vector
##'
##' @param x a vector
##' @return the midpoints, a vector of length length(x)-1
##' @export

midpts <- function(x){
    diffs <- diff(x)
    return(x[1:(length(x)-1)]+diffs/2)
}


##' getBbasis function
##'
##' A function returning the piecewise polynomial coefficients for a B-spline basis function i.e. the basis functions.
##'
##' @param x a vector of data
##' @param knots a vector of knots in ascending order. The first and last knots must be respectively the minimum and maximum of x.
##' @param degree the degree of the spline
##' @param force logical: skip check on knots? (not recommended!)
##' @return the knots and the piecewise polynomial coefficients for a B-spline basis function i.e. the basis functions.
##' @export

getBbasis <- function(x,knots,degree,force=FALSE){

    if(!all(sort(knots)==knots)){
        stop("Knots must be in ascending order")
    }
    if(!force){
        if(knots[1]!=min(x) | rev(knots)[1]!=max(x)){
            stop("First and last knots must be respectively the minimum and maximum of x")
        }
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



##' Bspline.construct function
##'
##' A function to construct a B-spline basis matrix for given data and basis coefficients. Used in evaluating the baseline hazard.
##'
##' @param x a vector, the data
##' @param basis an object created by the getBbasis function
##' @return a basis matrix
##' @export

Bspline.construct <- function(x,basis){
    idx <- as.numeric(cut(x,basis$knots,include.lowest=TRUE))
    xpows <- outer(x,0:(ncol(basis$poly[[1]])-1),"^")
    ans <- t(sapply(1:length(x),function(i){sapply(basis$poly,function(co){sum(co[idx[i],]*xpows[i,])})}))
    return(ans)
}



##' cumulativeBspline.construct function
##'
##' A function to construct the integral of a B-spline curve given data and basis coefficients. Used in evaluating the cumulative baseline hazard.
##'
##' @param x a vector, the data
##' @param basis an object created by the getBbasis function
##' @return an object that allows the integral of a given B-spline curve to be computed
##' @export

cumulativeBspline.construct <- function(x,basis){

    knots <- basis$knots

    coeffs <- basis$poly

    nb <- length(coeffs)
    np <- ncol(coeffs[[1]])

    addfun <- function(pars){ # if knot 1 is not zero, then need to add on integral from zero to first knot (add on rectangle -- continue first spline as a constant back to zero)
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


## cumulativeTPS.construct function
##
## A function to construct the integral of a B-spline curve given data and basis coefficients. Used in evaluating the cumulative baseline hazard.
##
## @param x a vector, the data
## @param basis an object created by the getBbasis function
## @return an object that allows the integral of a given B-spline curve to be computed
## @export
#
# FUNCTION IS UNFINISHED .NOT NECESSARY IF BOTH KNOTS INCLUDE 0 AS THE MINIMUM
#
# cumulativeTPS.construct <- function(x1,x2,basis1,basis2){
#
#     knots1 <- basis1$knots
#     knots2 <- basis2$knots
#
#     coeffs1 <- basis1$poly
#     coeffs2 <- basis2$poly
#
#     nb1 <- length(coeffs1)
#     np1 <- ncol(coeffs1[[1]])
#     nb2 <- length(coeffs2)
#     np2 <- ncol(coeffs2[[1]])
#
#     coeffs1 <- lapply(coeffs1,function(mat){return(t(apply(mat,1,function(x){x/(1:np)})))})
#     coeffs1 <- lapply(coeffs1,function(mat){return(cbind(0,mat))})
#     powers1 <- 0:np1
#     coeffs2 <- lapply(coeffs2,function(mat){return(t(apply(mat,1,function(x){x/(1:np)})))})
#     coeffs2 <- lapply(coeffs2,function(mat){return(cbind(0,mat))})
#     powers2 <- 0:np2
#
#     idx1 <- as.numeric(cut(x1,knots1,include.lowest=TRUE))
#     idx2 <- as.numeric(cut(x2,knots2,include.lowest=TRUE))
#
#     idxmax1 <- max(idx1,na.rm=TRUE)
#     idxmax2 <- max(idx2,na.rm=TRUE)
#
#     bm1 <- Bspline.construct(x=knots1,basis=basis1)
#     bm2 <- Bspline.construct(x=knots2,basis=basis2)
#
#     addfun <- function(mat){  # if knot 1 is not zero, then need to add on integral from zero to first knot (treat as constant to "left" of first knot)
#         #cof <- Reduce("+",mapply("*",basis$poly,pars,SIMPLIFY=FALSE))
#         cof <- matrix(0,np1,np2)
#         for(i in 1:np1){
#             for(j in 1:np2){
#                 cof[i,j] <- bm1[i,]%*%(m%*%bm2[j,])
#             }
#         }
#         if(knots1[1]==0 & knots2[1]==0){
#             return(0)
#         }
#         else{
#             return(knots[1]*sum(cof[1,]*knots[1]^(0:(np-1))))
#         }
#     }
#
#
#
#     ints1 <- t(sapply(coeffs1,function(mat){c(0,sapply(2:idxmax1,function(i){sum(mat[i-1,]*(knots1[i]^powers1-knots1[i-1]^powers1))}))}))
#     cumints2 <- t(apply(ints2,1,cumsum2))
#     ints2 <- t(sapply(coeffs2,function(mat){c(0,sapply(2:idxmax2,function(i){sum(mat[i-1,]*(knots2[i]^powers2-knots2[i-1]^powers2))}))}))
#     cumints2 <- t(apply(ints2,1,cumsum2))
#
#     integ <- c()
#     for (i in 1:length(idx)){
#         xpow <- x[i]^powers
#         kpow <- knots[idx[i]]^powers
#         integ <- rbind(integ,cumints[,idx[i]] + sapply(1:length(coeffs),function(j){sum(coeffs[[j]][idx[i],]*(xpow-kpow))}))
#     }
#
#     return(list(integral=integ,toadd=addfun))
# }
