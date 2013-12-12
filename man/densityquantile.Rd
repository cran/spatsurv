\name{densityquantile}
\alias{densityquantile}
\title{densityquantile function}
\usage{
densityquantile(obj, ...)
}
\arguments{
  \item{obj}{an object}

  \item{...}{additional arguments -- currently there are
  none, but this is for extensibility}
}
\value{
method densityquantile
}
\description{
Generic function for computing quantiles of the density
function for a given baseline hazard. This may not be
analytically tractable.
}
\seealso{
\link{densityquantile.basehazardspec},
\link{exponentialHaz}, \link{weibullHaz},
\link{gompertzHaz}, \link{makehamHaz}, \link{tpowHaz}
}

