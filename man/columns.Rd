\name{Get specific columns/rows fo a matrix}
\alias{columns}
\alias{rows}
\title{
Get specific columns/rows fo a matrix
}

\description{
Get specific columns/rows of a matrix.
}

\usage{
columns(x,indices)
rows(x,indices)
}

\arguments{
\item{x}{
A matrix with data.
}
\item{indices}{
An integer vector with the indices.
}
}


\value{
A matrix with the specific columns/rows of argumment indices.
}


\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}


\seealso{
\code{\link{rowMins}, \link{rowFalse}, \link{nth}, \link{colrange}, \link{colMedians}, \link{colVars}, \link{colSort}, \link{rowSort}, \link{rowTrue}
}
}

\examples{
x <- matrix(runif(100*100),100,100)
indices = sample(1:100,50)
all.equal(x[,indices],columns(x,indices))
all.equal(x[indices,],rows(x,indices))

x<-indices<-NULL
}

\keyword{ Extract columns/rows }