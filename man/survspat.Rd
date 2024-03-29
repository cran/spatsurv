% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survspat.R
\name{survspat}
\alias{survspat}
\title{survspat function}
\usage{
survspat(
  formula,
  data,
  dist,
  cov.model,
  mcmc.control,
  priors,
  shape = NULL,
  ids = list(shpid = NULL, dataid = NULL),
  control = inference.control(gridded = FALSE),
  boundingbox = NULL
)
}
\arguments{
\item{formula}{the model formula in a format compatible with the function flexsurvreg from the flexsurv package}

\item{data}{a SpatialPointsDataFrame object containing the survival data as one of the columns OR for polygonal data a data.frame, in which case, the argument shape must also be supplied}

\item{dist}{choice of distribution function for baseline hazard. Current options are: exponentialHaz, weibullHaz, gompertzHaz, makehamHaz, tpowHaz}

\item{cov.model}{an object of class covmodel, see ?covmodel ?ExponentialCovFct or ?SpikedExponentialCovFct}

\item{mcmc.control}{mcmc control parameters, see ?mcmcpars}

\item{priors}{an object of class Priors, see ?mcmcPriors}

\item{shape}{when data is a data.frame, this can be a SpatialPolygonsDataFrame, or a SpatialPointsDataFrame, used to model spatial variation at the small region level. The regions are the polygons, or they represent the (possibly weighted) centroids of the polygons.}

\item{ids}{named list entry shpid character string giving name of variable in shape to be matched to variable dataid in data. dataid is the second entry of the named list.}

\item{control}{additional control parameters, see ?inference.control}

\item{boundingbox}{optional bounding box over which to construct computational grid, supplied as an object on which the function 'bbox' returns the bounding box}
}
\value{
an object inheriting class 'mcmcspatsurv' for which there exist methods for printing, summarising and making inference from.
}
\description{
A function to run a Bayesian analysis on censored spatial survial data assuming a proportional hazards model using an adaptive Metropolis-adjusted
Langevin algorithm.
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
