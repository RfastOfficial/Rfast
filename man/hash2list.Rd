\name{Hash object to a list object}
\alias{hash2list}
\title{
Hash object to a list object
}

\description{
Hash object to a list object.
}

\usage{
hash2list(x, sorting = FALSE)
}

\arguments{
\item{x}{
A hash table with two parts, the keys (number(s) as string) and the key values (a single number). 
}
\item{sorting}{
This is if you you want the numbers in the keys sorted. The default value is FALSE.
}
}

\details{
For every key, there is a key value. This function creates a list and puts every pair of keys and value in a component of a list.
}

\value{
A list whose length is equal to the size of the hash table.
}

%\references{
%
%}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{  \link{hash.list}, \link{hash.find} }
}

\examples{
x=list("1 2 4 3"=2.56,"2.34 1.05"=2)
res<-hash2list(x)
res<-hash2list(x,TRUE)
}

\keyword{ Hash tables }
\keyword{ list }
