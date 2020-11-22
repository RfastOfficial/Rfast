\name{Cholesky decomposition of a square matrix}
\alias{cholesky}
\title{
Cholesky decomposition of a square matrix
}

\description{
Cholesky decomposition of a square matrix.
}

\usage{
cholesky(x,parallel = FALSE)
}

\arguments{
\item{x}{
A square positive definite matrix.
}
\item{parallel}{
A boolean value for parallel version.
}
}

\details{
The Cholesky decomposition of a square positive definite matrix is computed. The use of parallel is suggested for matrices with dimensions of 1000 or more.  
}

\value{
An upper triangular matrix.
}


\author{
Manos Papadakis

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm95@gmail.com>
}

\seealso{
\code{\link{is.symmetric}
}
}

\examples{
x = matrix(rnorm(1000 * 50), ncol = 50)
s = cov(x)
system.time(a1 <- cholesky(s))
system.time(a2 <- chol(s))
all.equal(a1[upper.tri(a1)], a2[upper.tri(a2)])
x <- NULL
s <- NULL
a1 <- NULL
a2 <- NULL
}

\keyword{ Cholesky decomposition }
