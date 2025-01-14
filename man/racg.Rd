\name{Angular central Gaussian random values simulation}
\alias{racg}
\title{
Angular central Gaussian random values simulation
}

\description{
Angular central Gaussian random values simulation.
}

\usage{
racg(n, sigma, seed = NULL)
}

\arguments{
\item{n}{
The sample size, a numerical value.
}
\item{sigma}{
The covariance matrix in \eqn{R^d}.
}
\item{seed}{
If you want the same to be generated again use a seed for the generator, an integer number.
}
}

\details{
The algorithm uses univariate normal random values and transforms them to multivariate via a spectral decomposition.
The vectors are then scaled to have unit length. 
}

\value{
A matrix with the simulated data.
}

\references{
Tyler D. E. (1987). Statistical analysis for the angular central Gaussian distribution on the sphere. Biometrika 74(3): 579-589. 
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> 
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{acg.mle}, \link{rmvnorm}, \link{rmvlaplace}, \link{rmvt}
}
}

\examples{
s <- cov( iris[, 1:4] )
#x <- racg(100, s)
#res<-acg.mle(x)  
#res<-vmf.mle(x)  ## the concentration parameter, kappa, is very low, close to zero, as expected.
}

\keyword{ Angular central Gaussian distribution }
\keyword{ random values simulation }
\keyword{ directional data }
