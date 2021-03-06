% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nonspatial_log_likelihood_and_derivatives.R
\name{maxlikparamPHsurv}
\alias{maxlikparamPHsurv}
\title{maxlikparamPHsurv function}
\usage{
maxlikparamPHsurv(surv, X, control)
}
\arguments{
\item{surv}{an object of class Surv}

\item{X}{the design matrix, containing covariate information}

\item{control}{a list containg various control parameters for the MCMC and post-processing routines}
}
\value{
initial estimates of the parameters
}
\description{
A function to get initial estimates of model parameters using maximum likelihood. Not intended for general purose use.
}
\references{
\enumerate{
    \item Benjamin M. Taylor and Barry S. Rowlingson (2017). spatsurv: An R Package for Bayesian Inference with Spatial Survival Models. Journal of Statistical Software, 77(4), 1-32, doi:10.18637/jss.v077.i04.
}
}
