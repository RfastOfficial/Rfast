\name{Column-row wise minima and maxima of two matrices}
\alias{colPmax}
\alias{colPmin}
\title{
Column-row wise minima and maxima of two matrices
}

\description{
Column-row wise minima and maxima of two matrices.
}

\usage{
colPmax(x, y)
colPmin(x, y)
}

\arguments{
\item{x}{
A numerical vector with numbers.
}
\item{y}{
A numerical vector with numbers.
}
}

\details{
The parallel minima or maxima are returned. This are the same as the base functions pmax and pmin. 
}

\value{
A numerical vector/matrix with numbers, whose length is equal to the length of the initital matrices containing the maximum or minimum between each pair. 
}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}

\seealso{
\code{ \link{Sort}, \link{colMins}, \link{colMaxs}, \link{colMedians}
}
}

\examples{
x <- matrix(rnorm(100),10,10)
y <- matrix(rnorm(100),10,10)
res<-colPmax(x, y)
res<-colPmin(x, y)
x<-y<-NULL
}

\keyword{ Pairs of vectors }
\keyword{ minimum }
\keyword{ maximum }
