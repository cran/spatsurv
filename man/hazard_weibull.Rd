\name{hazard_weibull}
\alias{hazard_weibull}
\title{hazard_weibull function}
\usage{
hazard_weibull(inputs)
}
\arguments{
  \item{inputs}{inputs for the function including the model
  matrix, frailties, fixed effects and the parameters of
  the baseline hazard derived from this model}
}
\value{
the hazard function for the individual
}
\description{
A function to compute the hazard function for an individual
where the baseline hazard comes from a Weibull survival
model
}

