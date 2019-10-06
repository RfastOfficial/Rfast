\name{Distance correlation}
\alias{dcor}
\alias{bcdcor}
\title{
Distance correlation
}

\description{
Distance correlation.
}

\usage{
dcor(x, y)
bcdcor(x, y)
}

\arguments{
\item{x}{
A numerical matrix.
}
\item{y}{
A numerical matrix.
}
}

\details{
The distance correlation or the bias corrected distance correlation of two matrices is calculated. The latter one 
is used for the hypothesis test that the distance correlation is zero (see \code{\link{dcor.ttest}}).
}

\value{
The value of the distance correlation of the bias corrected distance correlation.
}

\references{
G.J. Szekely, M.L. Rizzo and N. K. Bakirov (2007). Measuring and Testing Independence 
by Correlation of Distances. Annals of Statistics, 35(6):2769-2794.
}

\author{
Manos Papadakis

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Manos Papadakis <papadakm95@gmail.com>.
}

\seealso{
\code{ \link{dcov}, \link{dcor.ttest}, \link{edist}
}
}

\examples{
x <- as.matrix(iris[1:50, 1:4])
y <- as.matrix(iris[51:100, 1:4])
res<-dcor(x, y)
res<-bcdcor(x, y)

x<-y<-NULL
}

\keyword{ Distance correlation }
\keyword{ bias corrected }
