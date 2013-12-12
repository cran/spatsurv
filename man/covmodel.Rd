\name{covmodel}
\alias{covmodel}
\title{covmodel function}
\usage{
covmodel(model, pars)
}
\arguments{
  \item{model}{correlation type, see ?CovarianceFct}

  \item{pars}{vector of additional parameters for certain
  classes of covariance function (eg Matern), these must be
  supplied in the order given in ?CovarianceFct and are not
  estimated}
}
\value{
an object of class covmodel
}
\description{
A function to define the spatial covariance model, see also
?getcov and ?CovarianceFct
}

