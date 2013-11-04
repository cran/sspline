\name{smooth.sspline}
\alias{smooth.sspline}
\title{Smoothing Spline on the Sphere}
\description{
    It fits a smoothing splines on the sphere with the smoothing
    parameter chosen by the generalized cross validation (GCV) criteria
    or given by the user.
}
\usage{
  smooth.sspline(lon, lat, y, m = 2, smth = 0, lambda = 0)
}
\arguments{
  \item{lon}{ numeric vector, the longitudes }
  \item{lat}{ numeric vector, the latitudes }
  \item{y}{ numeric vector, the observations at (lon, lat) }
  \item{m}{ integer, order of smoothing, takes value from 1 to 10. Default to 2 }
  \item{smth}{ method for choosing the smoothing parameter: 0, gcv
    method; 1, user specified. Default to 0 }
  \item{lambda}{ used only when smth = 1. }
}
\details{
    It calls Fortran subroutine with the .Fortran interface.    
}
\value{
  A \code{smooth.sspline} object with the components
  \item{lon}{the original longitude}
  \item{lat}{the original latitude}
  \item{obs}{the original observation}
  \item{lambda}{the lambda that minimizes the gcv score}
  \item{gcv}{the corresponding gcv value at lambda}
  \item{varhat}{the estimated variance}
  \item{c}{the coefficient vector c for the estimated function}
  \item{d}{the coefficient d for the estimated function}
  \item{yhat}{the estimated (smoothed) observation}
  \item{call}{the call to \code{smooth.sspline}}
}
\references{Grace Wahba (1981), \emph{Spline Interpolation and Smoothing on the
    Sphere}, SIAM J. SCI. STAT. COMPUT.}
\author{ Xianhong Xie }
\note{The longitudes and latitudes are measured in degrees.}
\examples{
subdat <- WTdiff[sample(nrow(WTdiff), 200), 2:4]
attach(subdat)

smooth.sspline(lon, lat, avgd)

detach(subdat)
}
\keyword{smooth}


