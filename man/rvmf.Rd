\name{Simulation of random values from a von Mises-Fisher distribution}
\alias{rvmf}
\title{
Random values simulation from a von Mises-Fisher distribution
}

\description{
It generates random vectors following the von Mises-Fisher distribution. The data can be spherical or hyper-spherical.
}

\usage{
rvmf(n, mu, k, parallel = FALSE)
}

\arguments{
\item{n}{
The sample size.
}
\item{mu}{
The mean direction, a unit vector.
}
\item{k}{
The concentration parameter. If k = 0, random values from the spherical uniform will be drwan. 
Values from a multivariate normal distribution with zero mean vector and the identity matrix as 
the covariance matrix. Then each vector becomes a unit vector.
}
\item{parallel}{
Run the algorithm in parallel. The parallel algorithm needs memory.
}
}

\details{
It uses a rejection smapling as suggested by Andrew Wood (1994). 
}

\value{
A matrix with the simulated data.
}

\references{
Wood A. T. A. (1994). Simulation of the von Mises Fisher distribution. Communications in statistics-simulation 
and computation, 23(1): 157--164.

Dhillon I. S. & Sra S. (2003). Modeling data using directional distributions. Technical Report TR-03-06, 
Department of Computer Sciences, The University of Texas at Austin.
http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.75.4122&rep=rep1&type=pdf
}

\author{
Michail Tsagris and Manos Papadakis

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm85@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{vmf.mle}, \link{rvonmises}, \link{iag.mle}
}
}

\examples{
m <- rnorm(4)
m <- m/sqrt(sum(m^2))
x <- rvmf(1000, m, 25)
m
res<-vmf.mle(x)
}

\keyword{ Random values simulation }
\keyword{ von Mises-Fisher distribution }
