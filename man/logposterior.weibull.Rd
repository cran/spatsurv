\name{logposterior.weibull}
\alias{logposterior.weibull}
\title{logposterior.weibull function}
\usage{
logposterior.weibull(tm, delta, X, beta, omega, eta, gamma, priors, cov.model,
  u, control)
}
\arguments{
  \item{tm}{vector of observed times}

  \item{delta}{censoring indicator}

  \item{X}{design matrix}

  \item{beta}{beta vector at which to evaluate the
  posterior}

  \item{omega}{omegaa vector at which to evaluate the
  posterior}

  \item{eta}{eta vector at which to evaluate the posterior}

  \item{gamma}{the transformed latent Gaussian field at
  which to evaluate the posterior}

  \item{priors}{the priors, an object of class mcmcPriors}

  \item{cov.model}{an object of class covmodel, see
  ?covmodel}

  \item{u}{distance matrix}

  \item{control}{list of control parameters, see
  ?inference.control}
}
\value{
the log posterior

...
}
\description{
A function to evaluate the log-posterior
}

