\name{Minima and maxima of two vectors/matrices and Column-row wise minima and maxima of two matrices}
\alias{colPmax}
\alias{colPmin}
\alias{Pmax}
\alias{Pmin}
\alias{Pmin_Pmax}
\title{
Minima and maxima of two vectors/matrices and Column-row wise minima and maxima of two matrices
}

\description{
Minima and maxima of two vectors/matrices and Column-row wise minima and maxima of two matrices.
}

\usage{
colPmax(x, y)
colPmin(x, y)
Pmax(x, y,na.rm = FALSE)
Pmin(x, y,na.rm = FALSE)
Pmin_Pmax(x, y,na.rm = FALSE)
}

\arguments{
\item{x}{
A numerical vector or matrix with numbers.
}
\item{y}{
A numerical vector with numbers.
}
\item{na.rm}{
TRUE or FAlSE for remove NAs if exists.
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
