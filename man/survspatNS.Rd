% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survspatNS.R
\name{survspatNS}
\alias{survspatNS}
\title{survspatNS function}
\usage{
survspatNS(formula, data, dist, control = inference.control())
}
\arguments{
\item{formula}{the model formula in a format compatible with the function flexsurvreg from the flexsurv package}

\item{data}{a SpatialPointsDataFrame object containing the survival data as one of the columns}

\item{dist}{choice of distribution function for baseline hazard. Current options are: exponentialHaz, weibullHaz, gompertzHaz, makehamHaz, tpowHaz}

\item{control}{additional control parameters, see ?inference.control}
}
\value{
an object inheriting class 'mcmcspatsurv' for which there exist methods for printing, summarising and making inference from.
}
\description{
A function to perform maximun likelihood inference for non-spatial survival data.
}
\references{
\enumerate{
    \item Benjamin M. Taylor and Barry S. Rowlingson (2017). spatsurv: An R Package for Bayesian Inference with Spatial Survival Models. Journal of Statistical Software, 77(4), 1-32, doi:10.18637/jss.v077.i04.
}
}
\seealso{
\link{tpowHaz}, \link{exponentialHaz}, \link{gompertzHaz}, \link{makehamHaz}, \link{weibullHaz},
\link{covmodel}, \link{ExponentialCovFct}, \code{SpikedExponentialCovFct},
\link{mcmcpars}, \link{mcmcPriors}, \link{inference.control}
}
