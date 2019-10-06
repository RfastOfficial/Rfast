\name{Apply method to Positive and Negative number}
\alias{negative}
\alias{positive}
\alias{positive.negative}
\title{
Apply method to Positive and Negative number
}

\description{
Apply method to Positive and Negative number.
}

\usage{
negative(x,method = "min")
positive(x,method = "min")
positive.negative(x,method = "min")
}

\arguments{
\item{x}{
A numerical vector with data.
}
\item{method}{
Accept 3 values. "min", "max", "min.max".
}
}

\details{
These functions apply the chosen method to the chosen subset (negative, positive, or both) from the vector and return the result. 
}

\value{
negative: apply the chosen method to every negative number of the input vector.
positive: apply the chosen method to every positive number of the input vector.
positive.negative: apply the chosen method to every negative and positive number of the input vector.
}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}


\seealso{
\code{\link{nth}, \link{colnth}, \link{rownth},\link{sort_unique}, \link{Round}
}
}

\examples{
x <- rnorm(1000)

identical(negative(x,"min"), min(x<0))
identical(positive(x,"min"), min(x>0))
identical(positive.negative(x,"min"), c(min(x<0),min(x>0)))

x<-NULL
}

\keyword{ negative numbers }
\keyword{ positive numbers }
\keyword{ minimum or maximum of negative }
\keyword{ minimum or maximum of positive }