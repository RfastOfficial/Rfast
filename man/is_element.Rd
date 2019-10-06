\name{Find element}
\alias{is_element}
\title{
Find element
}

\description{
Search a value in an unordered vector.
}

\usage{
is_element(x, key)
}

\arguments{
\item{x}{
A vector or matrix with the data.
}
\item{key}{
A value to check if exists in the vector x.
}
}

\details{
Find if the key exists in the vector and return returns TRUE/FALSE if the value is been found. If the vector is unordered  it is fast but if the vector is ordered then use binary_search. The functions is written in C++ in order to be as fast as possible. 
}

\value{
TRUE/FALSE if the value is been found.
}

%\references{
%}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}


\seealso{
\code{ \link{binary_search} (buit-in R function)
}
}

\examples{
x <- rnorm(500)
key <- x[50]
b <- is_element(x, key) 
}

\keyword{ Find element }
\keyword{ Divide and Qonquer }
\keyword{ Linear time }
