\name{maxlikparamPHsurv}
\alias{maxlikparamPHsurv}
\title{maxlikparamPHsurv function}
\usage{
maxlikparamPHsurv(surv, X, control)
}
\arguments{
  \item{surv}{an object of class Surv}

  \item{X}{the design matrix, containing covariate
  information}

  \item{control}{a list containg various control parameters
  for the MCMC and post-processing routines}
}
\value{
initial estimates of the parameters
}
\description{
A function to get initial estimates of model parameters
using maximum likelihood. Not intended for general purose
use.
}

