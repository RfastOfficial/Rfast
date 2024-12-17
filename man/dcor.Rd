\name{Distance correlation}
\alias{dcor}
\title{
Distance correlation
}

\description{
Distance correlation.
}

\usage{
dcor(x, y, bc = FALSE)
}

\arguments{
\item{x}{
A numerical matrix.
}
\item{y}{
A numerical matrix.
}
\item{bc}{
If you want the bias-corrected distance correlation set this equal to TRUE.
}
}

\details{
The distance correlation or the bias corrected distance correlation of two matrices is calculated. The latter one 
is used for the hypothesis test that the distance correlation is zero (see \code{\link{dcor.ttest}}).
}

\value{
For the bias corrected distance correlation its value only. 
For the distance correlation a list including:
\item{dcov}{
The distance covariance.
}
\item{dvarX}{
The distance variance of x.
}
\item{dvarY}{
The distance variance of Y.
}
\item{dcor}{
The distance correlation.
}
}

\references{
G.J. Szekely, M.L. Rizzo and N. K. Bakirov (2007). Measuring and Testing Independence 
by Correlation of Distances. Annals of Statistics, 35(6): 2769--2794.

Szekely G. J. and Rizzo M. L. (2023). The Energy of Data and Distance Correlation. 
Chapman and Hall/CRC.
}

\author{
Manos Papadakis

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm95@gmail.com>.
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
