\name{Column and row-wise medians of a matrix or median of a vector.}
\alias{colMedians}
\alias{rowMedians}
\alias{Median}
\alias{med}
\title{
Column and row-wise medians of a matrix or median of a vector.
}

\description{
Column and row-wise medians of a matrix or median of a vector.
}

\usage{
colMedians(x,na.rm = FALSE, parallel = FALSE, cores = 0)
rowMedians(x,na.rm = FALSE, parallel = FALSE, cores = 0)
Median(x,na.rm=FALSE)
med(x,na.rm=FALSE)
}

\arguments{
\item{x}{
A vector, matrix or data.frame with the data.
}
\item{parallel}{
Do you want to do it in parallel in C++? TRUE or FALSE.
}
\item{na.rm}{
TRUE or FAlSE for remove NAs if exists.
}
\item{cores}{
Number of cores to use for parallelism. Valid only when argument parallel is set to TRUE. 
Default value is 0 and it means the maximum supported cores.
}
}

\details{
The functions is written in C++ in order to be as fast as possible.
}

\value{
A vector with the column medians.
}

\author{
R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}

\seealso{
\code{\link{Median}, \link{colVars}, \link{colMeans} (buit-in R function)
}
}

\examples{
x <- matrix( rnorm(100 * 100), ncol = 100 )
a <- apply(x, 2, median) 
b1 <- colMedians(x) 
all.equal(as.vector(a), b1)

x<-a<-b1<-NULL
}

\keyword{ Column-wise medians }
\keyword{ Row-wise medians }
\keyword{ medians }
