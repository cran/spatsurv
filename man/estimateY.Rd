% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calibratemcmc.R
\name{estimateY}
\alias{estimateY}
\title{estimateY function}
\usage{
estimateY(X, betahat, omegahat, surv, control)
}
\arguments{
\item{X}{the design matrix containing covariate information}

\item{betahat}{an estimate of beta}

\item{omegahat}{an estimate of omega}

\item{surv}{an object of class Surv}

\item{control}{a list containg various control parameters for the MCMC and post-processing routines}
}
\value{
an estimate of Y, to be used in calibrating the MCMC
}
\description{
A function to get an initial estimate of Y, to be used in calibrating the MCMC. Not for general use
}
