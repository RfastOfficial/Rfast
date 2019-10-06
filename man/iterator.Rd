\name{Iterator}
\alias{iterator}
\alias{print.iterator}
\alias{Elem<-}
\alias{Elem}
\alias{Elem<-.iterator}
\alias{Elem.iterator}
\alias{==.iterator}
\alias{!=.iterator}
\title{
Iterator
}

\description{
A way to traverse a list, data.frame, matrix or vector.
}

\usage{
iterator(x,method="ceil",type="vector",by=1)
\method{print}{iterator}(x,...)
\method{Elem}{iterator}(x) <- value
Elem(x)
Elem(x) <- value
\method{Elem}{iterator}(x)
\method{==}{iterator}(x,y)
\method{!=}{iterator}(x,y)
}

\arguments{
\item{x}{
A variable with any type, or iterator object.
}
\item{value}{
An value depending the method of the iterator.
}
\item{y}{
An iterator.
}
\item{method}{
Method of the iterator class. One of "ceil","col","row".
}
\item{type}{
One of "vector","matrix","data.frame","list".
}
\item{by}{
An integer value to iterate through element.
}
\item{...}{
Anything the user want.
}
}

\details{
iterator: is an object that helps a programmer to traverse the given object.

print.iterator: print an object of class iterator.

"Elem<-": access to element and change the value.

Elem: access to element.
}

\value{
An object of class "iterator". This object holds 4 fields:

copy: deep copy of iterator.
end: get iterator tha have access to points to the last element.
equals: equality of iterators
nextElem: move iterator to point to the next element using argument "by".
prevElem: move iterator to point to the previous element using argument "by".
}

\author{
R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}


\seealso{
\code{\link{colShuffle}, \link{colVars}, \link{colmeans}, \link{read.directory}
}
}

\examples{

y<-rnorm(100)
x<-iterator(y,method="ceil",type="vector",by=1)

s<-0
while(x != x$end()){
	s <- s + Elem(x)
	x$nextElem()
}

all.equal(s,sum(y))
}

\keyword{ iterator }
\keyword{ Environment }
