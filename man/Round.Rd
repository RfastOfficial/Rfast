\name{Round each element of a matrix/vector}
\alias{Round}
\title{
Round each element of a matrix/vector
}

\description{
Round each element of a matrix/vector.
}

\usage{
Round(x,digit=0,na.rm = FALSE)
}

\arguments{
\item{x}{
A numeric matrix/vector with data or NA. NOT integer values. 
}
\item{digit}{
An integer value for 0...N-1 where N is the number of the digits. By default is 0.
}
\item{na.rm}{
TRUE or FAlSE for remove NAs if exists.
}
}

\details{
Round is a very fast C++ implementation. Especially for large data. It handles NA.
}

\value{
A vector/matrix where each element is been rounded in the given digit.
}


\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}


\seealso{
\code{\link{Lchoose}, \link{Log}, \link{Choose}
}
}

\examples{
x <-matrix( rnorm( 50 * 10), ncol = 10 )
a <- Round(x,5)
b <- round(x,5)
all.equal(a,b) #true
x <-rnorm( 1000)
a <- Round(x,5)
b <- round(x,5)
all.equal(a,b) # true
}

\keyword{ Round vector/matrix }