\name{ Replicate columns/rows }
\alias{rep_col}
\alias{rep_row}
\title{
Replicate columns/rows
}

\description{
Replicate columns/rows.
}

\usage{
rep_col(x,n)
rep_row(x,n)
}

\arguments{
\item{x}{
A vector with data.
}
\item{n}{
Number of new columns/rows.
}
}


\value{
A matrix where each column/row is equal to "x".
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
x <- runif(10)
all.equal(rep_col(x,10),matrix(x,nrow=length(x),ncol=10))
all.equal(rep_row(x,10),matrix(x,ncol=length(x),nrow=10,byrow=TRUE))
}

\keyword{ Replicate in columns/rows }