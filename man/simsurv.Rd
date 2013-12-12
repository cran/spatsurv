\name{simsurv}
\alias{simsurv}
\title{simsurv function}
\usage{
simsurv(X = cbind(age = runif(100, 5, 50), sex = rbinom(100, 1, 0.5), cancer =
  rbinom(100, 1, 0.2)), beta = matrix(c(0.0296, 0.0261, 0.035), 3, 1),
  theta = 1, dist = "exp", coords = matrix(runif(2 * nrow(X)), nrow(X),
  2), sigmaphi = c(1, 0.1), cov.model = covmodel(model = "exponential", pars
  = NULL), mcmc.control = mcmcpars(nits = 1e+05, burn = 10000, thin = 90),
  savechains = TRUE)
}
\arguments{
  \item{X}{a matrix of covariate information}

  \item{beta}{the parameter effects}

  \item{theta}{parameter for the baseline hazard model (the
  rate for exponential data and the shape and scale for
  Weibull data)}

  \item{dist}{the distribution choice: exp or weibull at
  present}

  \item{coords}{matrix with 2 columns giving the
  coordinates at which to simulate data}

  \item{sigmaphi}{a vector of length 2: the parameters for
  the covariance function, sigma and phi in that order}

  \item{cov.model}{an object of class covmodel, see
  ?covmodel}

  \item{mcmc.control}{mcmc control paramters, see
  ?mcmcpars}

  \item{savechains}{save all chains? runs faster if set to
  FALSE, but then you'll be unable to conduct
  convergence/mixing diagnostics}
}
\value{
simulated survival times from the exponential model (the
last simulated value from the MCMC chains)
}
\description{
A function to simulate spatial parametric proportional
hazards model with baseline hazard derived from the
exponential or weibull model. The function works by
simulating candidate survival times using MCMC in parallel
for each individual based on each individual's covariates
and the common parameter effects, beta.
}
\seealso{
\link{covmodel}, \link{survspat}
}

