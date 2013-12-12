\name{plotsurv}
\alias{plotsurv}
\title{plotsurv function}
\usage{
plotsurv(spp, ss, maxcex = 1, background = NULL, eventpt = 19,
  eventcol = "red", censpt = "+", censcol = "black", xlim = NULL,
  ylim = NULL, ...)
}
\arguments{
  \item{spp}{A spatial points data frame}

  \item{ss}{A Surv object (with right-censoring)}

  \item{maxcex}{maximum size of dots default is equavalent
  to setting cex equal to 1}

  \item{background}{a background object to plot default is
  null, which gives a blamk background note that if
  non-null, the parameters xlim and ylim will be derived
  from this object.}

  \item{eventpt}{The type of point to illustrate events,
  default is 19 (see ?pch)}

  \item{eventcol}{the colour of events, default is black}

  \item{censpt}{The type of point to illustrate events,
  default is "+" (see ?pch)}

  \item{censcol}{the colour of censored observations,
  default is red}

  \item{xlim}{optional x-limits of plot, default is to
  choose this automatically}

  \item{ylim}{optional y-limits of plot, default is to
  choose this automatically}

  \item{...}{other arguments to pass to plot}
}
\value{
Plots the survival data non-censored observations appear as
dots and censored observations as crosses. The size of the
dot is proportional to the observed time.
}
\description{
A function to produce a 2-D plot of right censored spatial
survival data.
}

