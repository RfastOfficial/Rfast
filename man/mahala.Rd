\name{Mahalanobis distance}
\alias{mahala}
\title{Mahalanobis distance
}

\description{Mahalanobis distance.
}

\usage{
mahala(x, mu, sigma, ischol = FALSE)
}

\arguments{
\item{x}{
A matrix with the data, where rows denotes observations (vectors) and the columns contain the variables.
}
\item{mu}{
The mean vector.
}
\item{sigma}{
The covariance or any square symmetric matrix.
}
\item{ischol}{
A boolean variable set to true if the Cholesky decomposition of the covariance matrix is supplied in the argument \"sigma\".
}

}

%\details{
%For every value of \eqn{\alpha} the normal likelihood (see the refernece) is computed. At the end, the plot of the values is constructed.
%}

\value{
A vector with the Mahalanobis distances.
}

%\references{
%}

\author{
Matteo Fasiolo <matteo.fasiolo@gmail.com>, 

C++ and R implementation and documentation: Matteo Fasiolo <matteo.fasiolo@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{dista}, \link{colmeans}
}
}

\examples{
x <- matrix( rnorm(100 * 50), ncol = 50 )
m <- colmeans(x)
s <- cov(x)
a1 <- mahala(x, m, s) 
}

\keyword{ Mahalanobis distance }

