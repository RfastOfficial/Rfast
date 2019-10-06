\name{Hash - Pair function}
\alias{hash.list}
\title{
Hash - Pair function
}

\description{
Hash - Pair function.
}

\usage{
hash.list(key,x)
}

\arguments{
\item{key}{
The keys of the given values.
}
\item{x}{
The values.
}

}

\details{
This function pairs each item of of key and value make a unique hash table.
}

\value{
Returns the hash-list table.
}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{hash.find}
}
}

\examples{
x <- hash.list(letters,c(1:26))
x[["a"]]==1
}

\keyword{ Hash Function }
\keyword{ Pair Function}
