\name{Index of the columns of a data.frame which are a specific type}
\alias{which.is}
\title{
Index of the columns of a data.frame which are a specific type
}

\description{
Index of the columns of a data.frame which are a specific type.
}

\usage{
which.is(x,method="factor")
}

\arguments{
\item{x}{
A data.frame where some columns are expected to be factor variables.
}
\item{method}{
A character value about the type. One of, "numeric","factor","integer","logical".
}
}

\details{
The function is written in C++ and this is why it is very fast. 
}

\value{
A vector with the column indices which are factor variables. If there are no factor variables it will return an empty vector. 
}

\author{
Manos Papadakis <papadakm95@gmail.com>

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}

\seealso{
\code{ \link{nth}, \link{Match}
}
}

\examples{
res<-which.is(iris)
}

\keyword{ Factor variables}
\keyword{ Integer variables}
\keyword{ Numeric variables}
\keyword{ Logical variables}
\keyword{ Data Frame }
