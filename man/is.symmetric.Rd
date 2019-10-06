\name{Check whether a square matrix is symmetric}
\alias{is.symmetric}
\title{
Check whether a square matrix is symmetric
}

\description{
Check whether a square matrix is symmetric.
}

\usage{
is.symmetric(x)
}

\arguments{
\item{x}{
A square matrix with data. 
}
}

\details{
Instead of going through the whole matrix, the function will stop if the first disagreement is met. 
}

\value{
A boolean value, TRUE of FALSE.
}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}


\seealso{
\code{\link{cholesky}, \link{cora}, \link{cova}
}
}

\examples{
x <-matrix( rnorm( 100 * 400), ncol = 400 )
s1 <- cor(x)
is.symmetric(s1)
x <- x[1:100, ]
is.symmetric(x)

x<-s1<-NULL
}

\keyword{ Symmetric matrix }