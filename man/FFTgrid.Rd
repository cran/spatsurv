\name{FFTgrid}
\alias{FFTgrid}
\title{FFTgrid function}
\usage{
FFTgrid(spatialdata, cellwidth, ext)
}
\arguments{
  \item{spatialdata}{a SpatialPixelsDataFrame object}

  \item{cellwidth}{width of computational cells}

  \item{ext}{multiplying constant: the size of the extended
  grid: ext*M by ext*N}
}
\value{
a list
}
\description{
A function to generate an FFT grid and associated
quantities including cell dimensions, size of extended
grid, centroids,
}

