\name{Row-wise minimum and maximum}
\alias{rowMins}
\alias{rowMaxs}
\alias{rowMinsMaxs}
\title{
Row-wise minimum and maximum of a matrix
}

\description{
Row-wise minimum and maximum of a matrix.
}

\usage{
rowMins(x, value = FALSE)
rowMaxs(x, value = FALSE)
rowMinsMaxs(x)
}

\arguments{
\item{x}{
A numerical matrix with data.
}
\item{value}{
If the value is FALSE it returns the indices of the minimum/maximum, otherwise it returns the minimum and maximum values.
}
}


\value{
A vector with the relevant values.
}


\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}


\seealso{
\code{\link{colMins}, \link{colMaxs}, \link{nth}, \link{rowrange} \link{colMedians}, \link{colVars}, \link{colSort}, \link{rowSort}
}
}

\examples{
x <- matrix( rnorm(10 * 10), ncol = 10 )

s1 <- rowMins(x)
#s2 <- apply(x, 1, min)

s1 <- rowMaxs(x)
#s2 <- apply(x, 1, max)

#s1 <- c(apply(x, 1, min),apply(x, 1, max) )
s2 <- rowMinsMaxs(x)

x<-s1<-s2<-NULL
}

\keyword{ Row-wise minimum }
\keyword{ row-wise maximum }

