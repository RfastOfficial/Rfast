\name{Distance variance and covariance}
\alias{dvar}
\alias{dcov}
\title{
Distance variance and covariance
}

\description{
Distance variance and covariances.
}

\usage{
dvar(x)
dcov(x, y)
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
The distance variance of a matrix or the distance covariance of two matrices is calculated.
}

\value{
The distance covariance or distance variance.
}

\references{
G.J. Szekely, M.L. Rizzo and N. K. Bakirov (2007). Measuring and Testing Independence 
by Correlation of Distances. Annals of Statistics, 35(6):2769-2794.
}

\author{
Manos Papadakis

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{dcor}, \link{edist}
}
}

\examples{
x <- as.matrix(iris[1:50, 1:4])
y <- as.matrix(iris[51:100, 1:4])
res<-dcov(x, y)
res<-dvar(x)
}

\keyword{ Distance variance }
\keyword{ Distance covariance }
