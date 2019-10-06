\name{Sort - Integer Sort - Sort a vector coresponding to another}
\alias{Sort}
\alias{Sort.int}
\alias{sort_cor_vectors}
\title{
Sort - Integer Sort - Sort a vector coresponding to another 
}

\description{
Fast sorting a vector.
}

\usage{
Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
Sort.int(x)
sort_cor_vectors(x, base, stable = FALSE, descending = FALSE)
}

\arguments{
\item{x}{
A numerical/integer/character vector with data.
}
\item{base}{
A numerical/character vector to help sorting the x.
}
\item{descending}{
A boolean value (TRUE/FALSE) for sorting the vector in descending order. By default sorts the vector in ascending.
}
\item{partial}{
This argument has two usages. The first is an index number for sorting partial the vector.
The second is a vector with 2 values, start and end c(start,end). Gives you a vector where the elements between 
start and end will be sorted only. Not character vector.
}
\item{stable}{
A boolean value (TRUE/FALSE) for choosing a stable sort algorithm. Stable means that discriminates on the same elements. Not character vector.
}
\item{na.last}{
Accept 4 values. TRUE, FALSE, NA, NULL.

TRUE/FALSE: for put NAs last or first.

NA:  for remove NAs completely from vector.

NULL: by default. Leave it like that if there is no NA values.
}
}

\details{
This function uses the sorting algorithm from C++. The implementation is very fast and highly optimised. 
Especially for large data.
}

\value{
Sort and Sort.int: The sorted vector.

sort_cor_vectors: The first argument but sorted acording to the second.
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
system.time( s1 <- Sort(x) )
system.time( s2 <- sort(x) )
all.equal(s1,s2) #true  but not if many duplicates.

system.time( s1 <- Sort(x,partial=100) )
system.time( s2 <- sort(x,partial=100) )
all.equal(s1,s2) #true


system.time( s1 <- Sort(x,stable=TRUE) )
system.time( s2 <- sort(x) )
all.equal(s1,s2) #true

x <- as.character(x)
system.time( s1 <- Sort(x) )
system.time( s2 <- sort(x) )
all.equal(s1,s2) #true

y <- runif(1000)
b <- sort_cor_vectors(x,y)

x<-rpois(100,100)
all.equal(Sort.int(x),sort.int(x))

x<-y<-y<-s1<-s2<-NULL
}

\keyword{ Sort }
\keyword{ Stable Sort }
\keyword{ Sort 2 vectors }
