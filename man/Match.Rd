\name{ Match }
\alias{Match}
\title{
Match
}

\description{
Return the positions of its first argument that matches in its second.
}

\usage{
Match(x,key=NULL)
}

\arguments{
\item{x}{
A numeric vector.
}
\item{key}{
The value/vector for searching in vector x. For now let it NULL. \bold{dont't use it!}.
}
}

\details{
This function implements the R's \"match\" function.
This version basicaly calculates the match(x,sort(unique(x))) for now. Do not use the argument key!
}

\value{
Returns the position/positions of the given key/keys in the x vector.
}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{match}
}
}

\examples{
y <- rnorm(100)
a <- Match(y)
b <-50%in%y
all.equal(as.vector(a),as.vector(b))
}

\keyword{ Match Function }