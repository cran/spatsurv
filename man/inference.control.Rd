\name{inference.control}
\alias{inference.control}
\title{inference.control function}
\usage{
inference.control(gridded = FALSE, cellwidth = NULL, ext = 2)
}
\arguments{
  \item{gridded}{logical. Whether to perform compuation on
  a grid. Default is FALSE. Note in version 0.9-1, this is
  still in the testing phase.}

  \item{cellwidth}{the width of computational cells to use}

  \item{ext}{integer the number of times to extend the
  computational grid by in order to perform compuitation.
  The default is 2.}
}
\value{
returns parameters to be used in the function survspat
}
\description{
A function to control inferential settings.
}
\seealso{
\link{survspat}
}

