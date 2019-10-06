\name{Transpose of a matrix}
\alias{transpose}
\title{
Transpose of a matrix
}

\description{
Transpose of a matrix.
}

\usage{
transpose(x)
}

\arguments{
\item{x}{
A numerical \bold{square} matrix with data.
}
}


\value{
The transposed matrix.
}

\references{
Gilbert Strang (2006). Linear Algebra and its Applications (4th edition). 
}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}


\seealso{
\code{\link{nth}, \link{colMaxs}, \link{colMins}, \link{colrange}
}
}

\examples{
x <- matrix( rnorm(500 * 500), ncol = 500, nrow=500 )
system.time( transpose(x) )
system.time( t(x)  )

x<-NULL
}

\keyword{ Transpose }
\keyword{ matrix }

