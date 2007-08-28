\name{station}
\alias{station}
\title{Distribution of the Stations on the World}
\description{
    It gives a simple illumination on how the given (lon, lat)
    pairs distributes on the world.
}
\usage{
station(lon=NULL, lat=NULL, pch=24, col="blue", bg="red", ...)
}
\arguments{
  \item{lon}{ numeric, the longitudes }
  \item{lat}{ numeric, the latitudes }
  \item{pch}{ the plotting symbol }
  \item{col}{ color value or name, the color used to draw the symbol }
  \item{bg}{ color value or name, the color used to fill the sumbol }
  \item{\dots}{ other plotting parameters }
}
\details{
    It calls the map.world to draw a world map.
}
\value{
    NULL
}
\author{ Xianhong Xie <\email{xie@stat.wisc.edu}> }
\note{ The longitudes and latitudes are measured in degrees. }
\seealso{ \code{\link{map.world}} }
\examples{
data(WTdiff)
subdat <- WTdiff[sample(nrow(WTdiff), 200), 2:3]
attach(subdat)

station(lon, lat)

detach(subdat)
}
\keyword{hplot}

