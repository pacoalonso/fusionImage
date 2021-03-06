% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pan.R
\docType{data}
\name{pan}
\alias{pan}
\title{Test panchromatic image from Cuenca}
\format{Two objects of class \code{"raster"} (\code{mis} and \code{pan})}
\source{
\href{https://www.um.es}{UMU}
}
\usage{
data(image)
}
\description{
Test panchromatic image from Cuenca
}
\examples{
data(image)
pca_fusion(mis, pan, method="bilinear", bits=12)
}
\references{
Me lo pasó Fulgen
}
\keyword{datasets}
