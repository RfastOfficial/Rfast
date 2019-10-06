\name{Find the given value in a hash table}
\alias{hash.find}
\title{
Find the given value in a hash table
}

\description{
Find the given value in a hash table or list.
}

\usage{
hash.find(x,key)
}

\arguments{
\item{x}{
A hash table or list.
}

\item{key}{
The key for searching the table.
}
}

\details{
This function search the given key.
}

\value{
If the given key exists return its value else returns 0.
}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{hash.list}
}
}

\examples{
x <- hash.list(letters,c(1:26))
value <- hash.find(x,"a")
x[["a"]]==value
}

\keyword{ Hash Function }
\keyword{ Find Value}
