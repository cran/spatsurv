\name{weibull_ltar}
\alias{weibull_ltar}
\title{weibull_ltar function}
\usage{
weibull_ltar(t, XbetaplusY, expXbetaplusY, theta)
}
\arguments{
  \item{t}{vector of times}

  \item{XbetaplusY}{the matrix product X times beta where
  is the design mastrix and beta are the covariate effects}

  \item{expXbetaplusY}{the exponential of the matrix
  product X times beta where is the design mastrix and beta
  are the covariate effects}

  \item{theta}{vector of length 2: shape and scale
  parameters for the weibull model}
}
\value{
the log-likelihood
}
\description{
A function to evaluate the likelihood for the parametric
proportional hazards model with baseline hazard derived
from the weibull model.
}

