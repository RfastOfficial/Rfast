\name{Natural Logarithm each element of a matrix}
\alias{Log}
\title{
Natural Logarithm each element of a matrix
}

\description{
Natural Logarithm each element of a matrix.
}

\usage{
Log(x, na.rm = FALSE)
}

\arguments{
\item{x}{
A matrix with data. 
}
\item{na.rm}{
A boolean value (TRUE/FALSE) for removing NA. 
}
}

\details{
The argument must be a matrix. For vector the time was the same as R's "log" function so we did not add it.
}

\value{
A matrix where each element is the natural logarithm of the given argument.
}


\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}


\seealso{
\code{\link{Lbeta}, \link{Lchoose}, \link{Choose}
}
}

\examples{
x <-matrix( runif( 100 * 100), ncol = 100 )
a <- log(x) 
b <- Log(x) 
all.equal(a, b) # true

x<-a<-b<-NULL
}

\keyword{ Log matrix }