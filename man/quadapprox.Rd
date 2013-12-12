\name{quadapprox}
\alias{quadapprox}
\title{quadapprox function}
\usage{
quadapprox(fun, xseq, yseq, ...)
}
\arguments{
  \item{fun}{a function}

  \item{xseq}{sequence of x-values defining grid on which
  to compute the approximation}

  \item{yseq}{sequence of y-values defining grid on which
  to compute the approximation}

  \item{...}{other arguments to be passed to fun}
}
\value{
a 2 by 2 matrix containing the curvature at the maximum and
the (x,y) value at which the maximum occurs
}
\description{
A function to compute the second derivative of a function
using a quadratic approximation to the function on a grid
of points defined by xseq and yseq. Also returns the local
maximum.
}

