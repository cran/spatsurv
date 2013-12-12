\name{posteriorcov}
\alias{posteriorcov}
\title{posteriorcov function}
\usage{
posteriorcov(x, probs = c(0.025, 0.5, 0.975), rmax = NULL, n = 100,
  plot = TRUE, ...)
}
\arguments{
  \item{x}{an object of class mcmcspatsurv}

  \item{probs}{vector of probabilities to be fed to
  quantile function}

  \item{rmax}{maximum distance in space to compute this
  distance up to}

  \item{n}{the number of points at which to evaluate the
  posterior covariance.}

  \item{plot}{whether to plot the result}

  \item{...}{other arguments to be passed to matplot
  function}
}
\value{
produces a plot of the posterior spatial covariance
function.
}
\description{
A function to produce a plot of the posterior covariance
function with upper and lower quantiles.
}
\seealso{
\link{print.mcmcspatsurv}, \link{quantile.mcmcspatsurv},
\link{summary.mcmcspatsurv}, \link{vcov.mcmcspatsurv},
\link{frailtylag1}, \link{spatialpars}, \link{hazardpars},
\link{fixedpars}, \link{randompars}, \link{baselinehazard},
\link{predict.mcmcspatsurv}, \link{priorposterior},
\link{MCE}, \link{hazardexceedance}
}

