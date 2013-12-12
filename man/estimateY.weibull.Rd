\name{estimateY.weibull}
\alias{estimateY.weibull}
\title{estimateY.weibull function}
\usage{
estimateY.weibull(X, betahat, omegahat, tm, delta, u, covmodel)
}
\arguments{
  \item{X}{the design matrix}

  \item{betahat}{an estimate of beta}

  \item{omegahat}{an estimate of omega}

  \item{tm}{vector of observed times}

  \item{delta}{censoring indicator}

  \item{u}{vecor of distances bettween points}

  \item{covmodel}{an object of class 'covmodel', see
  ?covmodel}
}
\value{
an estimate of the latent Gaussian field
}
\description{
A function to estimate Y assuming a parametric proportional
hazards model with baseline hazard derived from the Weibull
distribution.
}

