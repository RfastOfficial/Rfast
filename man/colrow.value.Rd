\name{ Check if any column or row is fill with values }
\alias{colrow.value}
\title{
Check if any column or row is fill with values
}

\description{
Check if any column or row is fill with values.
}

\usage{
colrow.value(x,value=0)
}

\arguments{
\item{x}{
A vector with data.
}
\item{value}{
A value to check.
}
}

\details{
Check all the column if any has all its elements equal to argument value. If found, return "TRUE". Otherwise continues with rows. If columns and rows hasn't any value vector then return "FALSE". Even if it returns "FALSE" that doesn't mean the determinant can't be value. It might be but if check before and found any value vector then for sure the determinant it'll be value.
}

\value{
A boolean value, "TRUE" if any column OR row is all filled with value. "FALSE" otherwise.
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
x <- matrix(runif(10*10),10,10)
res<-colrow.value(x) 

x<-NULL
}

\keyword{ Column-Row wise checking }
\keyword{ Determinant }