\name{Column-wise minimum and maximum }
\alias{colMins}
\alias{colMaxs}
\alias{colMinsMaxs}
\title{
Column-wise minimum and maximum of a matrix
}

\description{
Column-wise minimum and maximum of a matrix.
}

\usage{
colMins(x, value = FALSE, parallel = FALSE)

colMaxs(x, value = FALSE, parallel = FALSE)

colMinsMaxs(x, parallel = FALSE)
}

\arguments{
\item{x}{
A numerical matrix or data.frame with data.
}
\item{value}{
If the value is FALSE it returns the indices of the minimum/maximum, otherwise it returns the minimum and maximum values.
}
\item{parallel}{
Do you want to do it in parallel in C++? TRUE or FALSE. The parallel will return the minimum/maximum value only. It will never return the indices.
}
}


\value{
A vector with the relevant values.
}

%\references{
%Shabalin, A.A. (2012). Matrix eQTL: Ultra fast eQTL analysis %via large matrix operations. Bioinformatics 28(10): 1353-1358.
%}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}

\seealso{
\code{\link{rowMins}, \link{rowMaxs}, \link{nth}, \link{colrange}, \link{colMedians}, \link{colVars}, \link{colSort}, \link{rowSort}
}
}

\examples{
x <- matrix( rnorm(100 * 200), ncol = 200 )

s1 <- colMins(x) 
s2 <- apply(x, 2, min) 

s1 <- colMaxs(x) 
s2 <- apply(x, 2, max) 

s1 <- colMinsMaxs(x)
s2 <- c(apply(x, 2, min), apply(x, 2, max)) 

x<-s1<-s2<-NULL
}

\keyword{ Column-wise minimum }
\keyword{ column-wise maximum }
\keyword{ column-wise minimum-maximum }

