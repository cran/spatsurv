\name{survspat}
\alias{survspat}
\title{survspat function}
\usage{
survspat(formula, data, dist, covmodel, mcmc.control, priors)
}
\arguments{
  \item{formula}{see ?flexsurvreg}

  \item{data}{a SpatialPointsDataFrame object}

  \item{dist}{choice of distribution function for baseline
  hazard options are: "exp"}

  \item{covmodel}{an object of class covmodel, see
  ?covmodel}

  \item{mcmc.control}{mcmc control parameters, see
  ?mcmcpars}

  \item{priors}{an object of class Priors, see ?mcmcPriors}
}
\value{
the mcmc output
}
\description{
A function to run a Bayesian analysis on right censored
survial data assuming a proportional hazards model with
baseline hazard derived from the exponential model.
}

