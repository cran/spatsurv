% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/postprocessing.R
\name{CSplot}
\alias{CSplot}
\title{CSplot function}
\usage{
CSplot(mod, plot = TRUE, bw = FALSE, ...)
}
\arguments{
\item{mod}{an object produced by the function survspat}

\item{plot}{whether to plot the result, default is TRUE}

\item{bw}{Logical. Plot in black/white/greyscale? Default is to produce a colour plot. Useful for producing plots for journals that do not accept colour plots.}

\item{...}{other arguments to pass to plot}
}
\value{
the x and y values used in the plot
}
\description{
A function to produce a diagnostic plot for model fit using the Cox-Snell residuals.
}
