\name{estimateY.exp}
\alias{estimateY.exp}
\title{estimateY.exp function}
\usage{
estimateY.exp(X, betahat, omegahat, tm, delta, u, covmodel)
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
hazards model with baseline hazard derived from the
exponential distribution.
}

