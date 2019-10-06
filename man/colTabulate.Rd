\name{Column and row-wise tabulate}
\alias{colTabulate}
\alias{rowTabulate}
\title{
Column and row-wise tabulate
}

\description{
Column and row-wise tabulate of a matrix.
}

\usage{
colTabulate(x, max_number = max(x))
rowTabulate(x, max_number = max(x))
}

\arguments{
\item{x}{
An integer matrix with the data. The numbers must start from 1, i.e. 1, 2, 3, 4,... No zeros are allowed. Anything else may cause a crash.
}
\item{max_number}{
The maximum value of vector x. If you know which is the max number use this argument for faster results or by default max(x). 
}
}

\details{
The functions is written in C++ in order to be as fast as possible.
}

\value{
A matrix where in each column the command "tabulate" has been performed. The number of rows of the returned matrix will be equal to the max_number if given. Otherwise, the functions
will find this number. 
}

%\references{
%}

\author{
R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}


\seealso{
\code{\link{colShuffle}, \link{colVars}, \link{colmeans} 
}
}

\examples{
x <- matrix( rbinom(100 * 100, 4, 0.5), ncol = 100 )
system.time( colTabulate(x) )
x <- t(x)
system.time( rowTabulate(x) )

x<-NULL
}

\keyword{ Column-wise tabulate }
\keyword{ Row-wise tabulate }
