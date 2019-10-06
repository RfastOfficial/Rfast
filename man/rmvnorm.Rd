\name{Multivariate normal and t random values simulation}
\alias{rmvnorm}
\alias{rmvt}
\title{
Multivariate normal and t random values simulation
}

\description{
Multivariate normal and t random values simulation.
}

\usage{
rmvnorm(n, mu, sigma)
rmvt(n, mu, sigma, v)
}

\arguments{
  \item{n}{
  The sample size, a numerical value.
}
  \item{mu}{
  The mean vector in \eqn{R^d}.
}
  \item{sigma}{
  The covariance matrix in \eqn{R^d}.
}
  \item{v}{
  The degrees of freedom.
}
}

\details{
The algorithm uses univariate normal random values and transforms them to multivariate via a spectral decomposition.
It is faster than the command "mvrnorm" available from MASS, and it allows for singular covariance matrices.
}

\value{
A matrix with the simulated data.
}

\references{
Aitchison J. (1986). The statistical analysis of compositional data. Chapman & Hall.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> 
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{racg}, \link{rmvlaplace}, \link{rmvt}
}
}

\examples{
x <- as.matrix(iris[, 1:4])
m <- colmeans(x)
s <- cov(x)
y <- rmvnorm(1000, m, s)
res<-colmeans(y)
res<-cov(y)
y <- NULL
}

\keyword{ multivariate normal distribution }
\keyword{ random values simulation }
