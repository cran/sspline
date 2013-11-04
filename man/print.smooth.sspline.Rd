\name{print.smooth.sspline}
\alias{print.smooth.sspline}
\alias{summary.smooth.sspline}
\title{Display a Smooth.sspline Object}
\description{
    The print and summary methods for \code{smooth.sspline}
    object.
}
\usage{
  \method{print}{smooth.sspline}(x, ...)
  \method{summary}{smooth.sspline}(object, ...)
}
\arguments{
    \item{ x }{ \code{smooth.sspline} objects }
    \item{ object }{ \code{smooth.sspline} objects }
    \item{\dots}{ other parameters, not used }
}
\value{
    For \code{print.smooth.sspline}, a \code{smooth.sspline} object; for
    \code{summary.smooth.sspline}, NULL.
  }
\author{ Xianhong Xie }
\examples{
subdat <- WT6367[sample(nrow(WT6367), 200), 2:4]
attach(subdat)

splobj <- smooth.sspline(lon, lat, avgt)

print(splobj)
summary(splobj)

detach(subdat)
}
\keyword{smooth}
\keyword{methods}


