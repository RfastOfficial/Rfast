\name{Variance of a vector}
\alias{Var}
\title{
Variance (and standard deviation)  of a vector
}

\description{
Variance (and standard deviation)  of a vector.
}

\usage{
Var(x, std = FALSE,na.rm = FALSE)
}

\arguments{
\item{x}{
A vector with data. 
}
\item{std}{
If you want the standard deviation set this to TRUE, otherwise leave it FALSE.
}
\item{na.rm}{
TRUE or FAlSE for remove NAs if exists.
}
}

\details{
This is a faster calculation of the usual variance of a matrix. 
}

\value{
The variance of the vector.
}


\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Manos Papadakis <papadakm95@gmail.com>.
}


\seealso{
\code{\link{colVars}, \link{cova}
}
}

\examples{
x <- rnorm(100)
a1 <- Var(x)
a2 <- var(x)
x<-NULL
}

\keyword{ Variance }
