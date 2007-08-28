\name{predict.smooth.sspline}
\alias{predict.smooth.sspline}
\title{Spherical Smoothing Spline Prediction}
\description{
    Make prediction on the sphere using the information
    got from a \code{smooth.sspline} object.
}
\usage{
predict.smooth.sspline(object, lon, lat, grid=FALSE, ...)
}
\arguments{
    \item{object}{ a \code{smooth.sspline} object }
    \item{lon}{ the longitudes on which the prediction is to be made }
    \item{lat}{ the latitudes on which the prediction is to be made }
    \item{grid}{ whether the prediction is on a grid }
    \item{\dots}{ other parameters, not used }
}
\details{
    It calls Fortran subroutine with the .Fortran interface.
}
\value{
    If grid = TRUE, return a matrix with dimension (length(lon),
    length(lat)); otherwise, return a vector of length = length(lon).    
}
\references{Grace Wahba (1981), \emph{Spline Interpolation and Smoothing
	on the Sphere}, SIAM J. SCI. STAT. COMPUT. }
\author{ Xianhong Xie }
\note{ The longitudes and latitudes are measured in degrees. }
\seealso{\code{\link{smooth.sspline}}}

\examples{
data(WT9397)

subdat <- WT9397[sample(nrow(WT9397), 200), 2:4]
attach(subdat)

splobj <- smooth.sspline(lon, lat, avgt)

predict(splobj, lon=seq(-180,180,len=50), lat=seq(-90,90,len=25), grid=TRUE)

detach(subdat)
}
\keyword{smooth}
\keyword{methods}

