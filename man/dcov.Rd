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
A numerical matrix or a vector.
}
\item{y}{
A numerical matrix or a vector.
}
}

\details{
The distance variance of a matrix/vector or the distance covariance of two matrices is calculated.
For the distance variance of a vector we use the fast method of Huo and Szekely (2016).
}

\value{
The distance covariance or distance variance.
}

\references{
Szekely G.J., Rizzo M.L. and Bakirov N.K.(2007). Measuring and Testing Independence 
by Correlation of Distances. Annals of Statistics, 35(6):2769-2794.

Huo X. and Szekely G. J. (2016). Fast computing for distance covariance. 
Technometrics, 58(4): 435-447.
}

\author{
Manos Papadakis

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm95@gmail.com>.
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
res <- dcov(x, y)
res <- dvar(x[, 1])
}

\keyword{ Distance variance }
\keyword{ Distance covariance }
