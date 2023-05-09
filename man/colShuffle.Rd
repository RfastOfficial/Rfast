\name{Column and row-wise Shuffle}
\alias{colShuffle}
\alias{rowShuffle}
\title{
Column and row-wise Shuffle
}

\description{
Column and row-wise shuffle of a matrix.
}

\usage{
colShuffle(x)
rowShuffle(x)
}

\arguments{
\item{x}{
A matrix or data.frame with the data.
}
}

\details{
The functions is written in C++ in order to be as fast as possible.
}

\value{
A vector with the column/row Shuffle.
}

%\references{
%Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for compositional data. In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain.
%}

\author{
R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}


\seealso{
\code{\link{Median}, \link{colVars}, \link{colMeans} (buit-in R function)
}
}

\examples{
x <- matrix( rnorm(100 * 100), ncol = 100 )
system.time( colShuffle(x) )
system.time( rowShuffle(x) )

x<-NULL
}

\keyword{ Column-wise Shuffle }
\keyword{ Row-wise Shuffle }
