\name{ Sort and unique numbers }
\alias{sort_unique}
\alias{sort_unique.length}
\title{
Sort and unique
}

\description{
Sort and unique numbers.
}

\usage{
sort_unique(x)
sort_unique.length(x)
}

\arguments{
\item{x}{
A numeric vector. 
}
}

\details{
The "sort_unique" function implements R's "unique" function using C++'s function but also sort the result. The "sort_unique.length" returns the length of the unique numbers only for \bold{itegers}. 
}

\value{
Returns the discrete values but sorted or their length (depending on the function you do).
}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{colSort}, \link{rowSort}, \link{sort_cor_vectors}
}
}

\examples{
y <- rnorm(100)
a <- sort_unique(y)
b <- sort.int(unique(y))
all.equal(as.vector(a),as.vector(b))
x <- rpois(1000,10)
sort_unique.length(x)
length(sort_unique(x))

x<-a<-b<-NULL
}

\keyword{ Sort function }
\keyword{ unique numbers }