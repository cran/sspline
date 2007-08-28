\name{WT6367}
\alias{WT6367}
\non_function{}
\title{World Average Winter Temperature from 1963-1967}
\description{
    The \code{WT6367} data frame has 1391 rows and 4 columns.
    It contains the average temperature from 1963 to 1967 for
    those stations having non-missing observations on the
    winter (Dec-Feb) for ten years (1963-1967 and 1993-1997).
}
\format{
    This data frame contains the following columns:
    \describe{
	\item{recid}{
	    a numeric vector containing the coded information of the
	    stations (length 11). The first three digits represent the
	    country code; the next five digits, the station number; the
	    last three digits, whether a station is a WMO staion or close
	    to one.
	}
	\item{lon}{
	    a numeric vector containing the longitudes (in degrees) of
	    the stations.
	}
	\item{lat}{
	    a numeric vector containing the latitudes (in degrees) of
	    the stations.
	}
	\item{avgt}{
	    a numeric vector containing the average temperatures for
	    the stations (rounded to the second decimal point). 
	}
    }
}
\source{
    The Global Historical Climatology Network (GHCN) \cr
    \url{http://www.ncdc.noaa.gov/cgi-bin/res40.pl?page=ghcn.html}
}
\examples{
data(WT6367)

## Fit a smoothing spherical spline with part of the data
subdat <- WT6367[sample(nrow(WT6367), 200), 2:4]
attach(subdat)

smooth.sspline(lon, lat, avgt)

detach(subdat)
}
\keyword{datasets}


