\name{Vector allocation in a symmetric matrix}
\alias{squareform}
\title{
Vector allocation in a symmetric matrix
}

\description{
Vector allocation in a symmetric matrix.
}

\usage{
squareform(x)
}

\arguments{
\item{x}{
An numverical vector whose size must be the one that matches the dimensions of the final matrix. See examples.
}
}

\details{
The functions is written in C++ in order to be as fast as possible.
}

\value{
A symmetric matrix. The vector is allocated in the upper and in the lower part of the matrix. The diagonal is filled
with zeros.
}


\author{
R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}


\seealso{
\code{\link{colShuffle}, \link{colVars}, \link{colmeans} 
}
}

\examples{
x <- rnorm(1)
res<-squareform(x)  ## OK
x <- rnorm(3)
res<-squareform(x)  ## OK
x <- rnorm(4)
res<-squareform(x)  ## not OK
}

\keyword{ symmetric matrix }
