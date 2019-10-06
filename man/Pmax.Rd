\name{Minima and maxima of two vectors/matrices}
\alias{Pmax}
\alias{Pmin}
\alias{Pmin_Pmax}
\title{
Minima and maxima of two vectors/matrices
}

\description{
Minima and maxima of two vectors/matrices.
}

\usage{
Pmax(x, y,na.rm = FALSE)
Pmin(x, y,na.rm = FALSE)
Pmin_Pmax(x, y,na.rm = FALSE)
}

\arguments{
\item{x}{
A numerical vector with numbers.
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
A numerical vector/matrix with numbers, whose length is equal to the length of the initial vectors/matrices containing the maximum or minimum between each pair. 
}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}

\seealso{
\code{ \link{colSort}, \link{rowSort}, \link{Sort}, \link{colMins}
}
}

\examples{
x <- rnorm(10)
y <- rnorm(10)
res<-Pmax(x, y)
a<-pmax(x, y)
res<-Pmin(x, y)
b<-pmin(x, y)
res<-Pmin_Pmax(x,y) == c(a,b)
a<-b<-x<-y<-NULL
}

\keyword{ Pairs of vectors }
\keyword{ minimum }
\keyword{ maximum }
