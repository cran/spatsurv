##' spatsurv
##'
##' An R package for spatially correlated parametric proportional hazards survial analysis.
##'
##' \packageDESCRIPTION{spatsurv}
##' \packageIndices{spatsurv}
##'
##'
##' section{Dependencies}{
##' The package \code{spatsurv} depends upon some other important contributions to CRAN in order to operate; their uses here are indicated:\cr\cr
##'     survival, sp, spatstat, raster, iterators, RandomFields, fields, rgl, Matrix, stringr, RColorBrewer, geostatsp, rgeos.
##' }
##'
##' section{Citation}{
##' To cite use of \code{spatsurv}, the user may refer to the following work:\cr\cr
##' Benjamin M. Taylor and Barry S. Rowlingson (2017).\cr
##' spatsurv: An R Package for Bayesian Inference with Spatial Survival Models.\cr
##' Journal of Statistical Software, 77(4), 1-32, doi:10.18637/jss.v077.i04.
##' }
##'
##' references{
##' X
##' }
##'
##' @docType package
##' @name spatsurv-package
##' @author Benjamin Taylor, Health and Medicine, Lancaster University,
##'  Barry Rowlingson, Health and Medicine, Lancaster University
##' @keywords package
##'
##'



## @importFrom OpenStreetMap openmap
##' @importFrom RColorBrewer brewer.pal
##' @importFrom stringr str_count str_detect
##' @importFrom Matrix nearPD Matrix sparseMatrix bdiag
## @importFrom rgl abclines3d aspect3d axes3d planes3d points3d segments3d text3d title3d
##' @importFrom fields image.plot
##' @importFrom RandomFields CovarianceFct
##' @importFrom rgeos gBuffer
##' @importFrom iterators icount iter nextElem
##' @importFrom sp coordinates<- bbox proj4string<- proj4string SpatialPixelsDataFrame SpatialGridDataFrame Polygon Polygons SpatialPolygons coordinates CRS geometry GridTopology over proj4string SpatialGrid SpatialPixels SpatialPoints SpatialPolygonsDataFrame split spTransform
#' @importFrom spatstat.core rpoint
#' @importFrom spatstat.geom inside.owin progressreport
##' @importFrom survival Surv survfit
##' @importFrom raster crop brick raster
##' @import lubridate


## @import stats
##' @import graphics
## @import methods
## @import utils
## @import grDevices

##' @importFrom stats vcov as.formula acf coefficients deriv dexp dist dnorm end fft fitted formula Gamma integrate knots lm model.matrix optim optimise poly quantile rbinom rexp rnorm runif sd start update var residuals cov
##' @importFrom graphics boxplot polygon hist legend lines matplot par plot points title abline
##' @importFrom methods as
##' @importFrom utils txtProgressBar setTxtProgressBar browseURL flush.console
##' @importFrom grDevices adjustcolor





`spatsurv` = NA
