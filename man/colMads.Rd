\name{Column and rows-wise mean absolute deviations}
\alias{colMads}
\alias{rowMads}
\title{
Column and row-wise mean absolute deviations
}

\description{
Column and row-wise mean absolute deviations.
}


\usage{
colMads(x,method = "median",na.rm=FALSE,parallel = FALSE)
rowMads(x,method = "median",na.rm=FALSE,parallel = FALSE)
}

\arguments{
\item{x}{
A matrix or data.frame with the data.
}
\item{method}{
A character vector with values "median", for median absolute deviation or "mean", for mean absolute deviation. 
}
\item{na.rm}{
A logical value TRUE/FALSE to remove NAs.
}
\item{parallel}{
A boolean value for parallel version.
}
}

\details{
The functions is written in C++ in order to be as fast as possible.
}

\value{
A vector with the column-wise mean absolute deviations.
}

%\references{
%Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for compositional data. In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain.
%}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm95@gmail.com>.
}


\seealso{
\code{\link{colMedians}, \link{rowMedians}, \link{colVars}, \link{colmeans}, \link{colMeans} (buit-in R function)
}
}

\examples{
x <- matrix( rnorm(100 * 100), ncol = 100 )
system.time( a <- colMads(x) )

x<-NULL
}

\keyword{ Column-wise median absolute deviations }