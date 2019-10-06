\name{Binary search algorithm}
\alias{binary_search}
\title{
Binary search algorithm
}

\description{
Search a value in an ordered vector. 
}

\usage{
binary_search(x, v, index=FALSE)
}

\arguments{
\item{x}{
A vector with the data.
}
\item{v}{
A value to check if exists in the vector x.
}
\item{index}{
A boolean value for choose to return the position inside the vector.
}
}

\details{
The functions is written in C++ in order to be as fast as possible.
}

\value{
Search if the v exists in x. Then returns TRUE/FALSE if the value is been found.
}

%\references{
%Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for compositional data. In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain.
%}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}


\seealso{
\code{\link{is_element} 
}
}

\examples{
x <- sort(rnorm(1000))
v <- x[50]
b <- binary_search(x,v) 
b1 <- binary_search(x,v,TRUE) 
}

\keyword{ Binary search Algorithm }
\keyword{ Divide and Qonquer }

