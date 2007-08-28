\name{plot.smooth.sspline}
\alias{plot.smooth.sspline}
\title{Plot a Smooth.sspline Object}
\description{
    Plot a smoothing spherical spline using color to represent
    the function value.
}
\usage{
plot.smooth.sspline(x, lon, lat, main="", xlab="Longitude",
    ylab="Latitude", key.title="Temp\n(deg)", ...)
}
\arguments{
  \item{x}{ a \code{smooth.sspline} object }
  \item{lon}{ the longitudes on which the function values will be calculated }
  \item{lat}{ the latitudes on which the function values will be calculated }
  \item{main}{ the main title of the plot }
  \item{xlab}{ the x-axis label of the main plot }
  \item{ylab}{ the y-axis label of the main plot }
  \item{key.title}{ the title for the colored key }
  \item{\dots}{ other plotting parameters, such as \emph{lwd}, \emph{asp}, and \ldots }
}
\details{
    It calls \code{predict.smooth.sspline} and \code{filled.contour}.
}
\value{
    NULL
}
\author{ Xianhong Xie }
\note{ The longitudes and latitudes are measured in degrees. }
\seealso{ \code{\link{predict.smooth.sspline}} }

\examples{
data(WTdiff)

subdat <- WTdiff[sample(nrow(WTdiff), 200), 2:4]
attach(subdat)

splobj <- smooth.sspline(lon, lat, avgd)

plot(splobj, lon=seq(-180, 180, len=50), lat=seq(-90, 90, len=25),
     main="World Average Temperature Change")

detach(subdat)
}
\keyword{smooth}
\keyword{methods}
\keyword{hplot}
