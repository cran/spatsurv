\name{exp_ltar}
\alias{exp_ltar}
\title{exp_ltar function}
\usage{
exp_ltar(t, XbetaplusY, expXbetaplusY, theta)
}
\arguments{
  \item{t}{vector of times}

  \item{XbetaplusY}{the matrix product X times beta where
  is the design mastrix and beta are the covariate effects}

  \item{expXbetaplusY}{the exponential of the matrix
  product X times beta where is the design mastrix and beta
  are the covariate effects}

  \item{theta}{rate parameter for exponential model}
}
\value{
the log-likelihood
}
\description{
A function to evaluate the likelihood for the parametric
proportional hazards model with baseline hazard derived
from the exponential model.
}

