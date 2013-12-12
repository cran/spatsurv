\name{transformestimates.weibull}
\alias{transformestimates.weibull}
\title{transformestimates.weibull function}
\usage{
transformestimates.weibull(x)
}
\arguments{
  \item{x}{a vector of paramters}
}
\value{
the transformed parameters. For the weibull model, this
transforms 'shape' 'scale' (see ?dweibull) to 'alpha' and
'lambda' for the MCMC
}
\description{
A function to transform estimates of the parameters of the
weibull baseline hazard function, so they are commensurate
with R's inbuilt density functions.
}

