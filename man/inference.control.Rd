\name{inference.control}
\alias{inference.control}
\title{inference.control function}
\usage{
inference.control(gridded = FALSE, cellwidth = NULL, ext = 2,
  MLinits = NULL, plotcal = FALSE)
}
\arguments{
  \item{gridded}{logical. Whether to perform compuation on
  a grid. Default is FALSE.}

  \item{cellwidth}{the width of computational cells to use}

  \item{ext}{integer the number of times to extend the
  computational grid by in order to perform compuitation.
  The default is 2.}

  \item{MLinits}{optional initial values for the
  non-spatial maximum likelihood routine used to initialise
  the MCMC, a vector of length equal to the number of
  parameters of the baseline hazard}

  \item{plotcal}{logical, whether to produce plots of the
  MCMC calibration process, this is a technical option and
  should onyl be set to TRUE if poor mixing is evident (the
  printed h is low), then it is also useful to use a
  graphics device with multiple plotting windows.}
}
\value{
returns parameters to be used in the function survspat
}
\description{
A function to control inferential settings. This function
is used to set parameters for more advanced use of
spatsurv.
}
\seealso{
\link{survspat}
}

