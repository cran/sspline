\name{map.world}
\alias{map.world}
\title{World Map}
\description{
    Sketch the continental boundary to give a rough idea of
    the position on the world.
}
\usage{
map.world(add=FALSE, main="", ...)
}
\arguments{
    \item{add}{ a logical scalar, if TRUE, add a map to the
	existing plot; otherwise, plot a new world map }
    \item{main}{ a character vector, the main title of the plot }
    \item{\dots}{ other parameters needed to pass to the lines
	function }
}
\value{
    NULL
}
\references{ S Archive under \url{http://lib.stat.cmu.edu} }
\author{ Original in S by Steve Wofsy <\email{scw@io.harward.edu}>, \cr
    ported to R by Xianhong Xie. }
\examples{
map.world(main = "The World Map")
}
\keyword{aplot}
\keyword{hplot}
