\name{Spherical and hyperspherical median}
\alias{mediandir}
\title{
Fast calculation of the spherical and hyperspherical median
}
\description{
It calculates, very fast, the (hyper-)spherical median of a sample.
}
\usage{
mediandir(x)
}
\arguments{
\item{x}{
The data, a numeric matrix with unit vectors.
}
}
\details{
The "mediandir" employes a fixed poit iterative algorithm stemming from the first derivative (Cabrera and Watson, 1990)
to find the median direction as described in Fisher (1985) and Fisher, Lewis and Embleton (1987).
}
\value{
The median direction.
}
\references{
Fisher N. I. (1985). Spherical medians. Journal of the Royal Statistical Society. Series B, 47(2): 342-348.

Fisher N. I., Lewis T. and Embleton B. J. (1987). Statistical analysis of spherical data. Cambridge university press.

Cabrera J. and Watson G. S. (1990). On a spherical median related distribution. Communications in Statistics-Theory and Methods, 19(6): 1973-1986.
}
\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> 
}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{ \link{vmf.mle} }
}
\examples{
m <- rnorm(3)
m <- m / sqrt( sum(m^2) )
#x <- rvmf(100, m, 10)
#res<-mediandir(x)
x <- NULL
}
\keyword{ Median direction }
