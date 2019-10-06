\name{Column-wise differences}
\alias{coldiffs}
\title{
Column-wise differences
}

\description{
Column-wise differences.
}

\usage{
coldiffs(x)
}

\arguments{
\item{x}{
A matrix with numbers.
}
}

\details{
This function simply does this function x[, -1] - x[, -k], where k is the last column of the matrix x. But it does it a lot faster. 
That is, 2nd column - 1st column, 3rd column - 2nd column, and so on.
}

\value{
A matrix with one column less containing the differences between the successive columns. 
}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}

\seealso{
\code{ \link{Dist}, \link{dista}, \link{colmeans}
}
}

\examples{
x <- matrix( rnorm(50 * 10), ncol = 10 )
res<-coldiffs(x)

x<-NULL
}

\keyword{ Differences }
\keyword{ matrix }
