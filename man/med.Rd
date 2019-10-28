\name{Median of a vector}
\alias{med}
\alias{Median}
\title{
Median of a vector
}

\description{
Median of a vector.
}

\usage{
med(x,na.rm=FALSE)
Median(x,na.rm=FALSE)
}

\arguments{
\item{x}{
A numerical vector.
}
\item{na.rm}{
TRUE or FAlSE for remove NAs if exists.
}
}

\details{
The function is written in C++ and this is why it is very fast.  
}

\value{
The median of the vector of a numbers.
}

\author{
Manos Papadakis <papadakm95@gmail.com>

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Manos Papadakis <papadakm95@gmail.com>.
}

\seealso{
\code{\link{nth}, \link{colMedians}
}
}

\examples{
x <- rnorm(1000)
a1 <- Median(x) 
a2 <- median(x)
}

\keyword{ nth elements }
\keyword{ sorting } 