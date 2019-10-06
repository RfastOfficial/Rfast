\name{Multivariate Laplace random values simulation}
\alias{rmvlaplace}
\title{
Multivariate Laplace random values simulation
}

\description{
Multivariate Laplace random values simulation.
}

\usage{
rmvlaplace(n, lam, mu, G) 
}

\arguments{
\item{n}{
The sample size, a numerical value.
}
\item{lam}{
The the parameter of the exponential distribution, a positive number.
}
\item{mu}{
The mean vector.
}
\item{G}{
A \eqn{d \times d} covariance matrix with determinant 1. 
}
}

\details{
The algorithm uses univariate normal random values and transforms them to multivariate via a spectral decomposition.
}

\value{
A matrix with the simulated data.
}

\references{
Eltoft T., Kim T., and Lee T.W. (2006). On the multivariate laplace distribution. Signal Processing Letters, IEEE, 13(5):300-303.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> 
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{rmvnorm}, \link{racg}, \link{rmvt}
}
}

\examples{
m <- colmeans( as.matrix( iris[, 1:4] ) )
s <- cov(iris[,1:4])
s <- s / det(s)^0.25
lam <- 3
x <- rmvlaplace(100, lam, m, s)
}

\keyword{ multivariate Laplace distribution }
\keyword{ random values simulation }