\name{ Represantation of Stack }
\alias{ Stack }
\title{
Represantation of Stack
}

\description{
Represantation of Stack.
}

\usage{
Stack(x,type=NULL)
}

\arguments{
\item{x}{
Any type that could be convert to vector or an integer value.
}
\item{type}{
A type for the Stack, "integer", "numeric" or any other that accepts one argument.
}
}

\details{
Stack is an abstract data type - data structure based on the principle of last in first out. 
To access the 3 fields, use operator "$".
}

\value{
An object of class "Stack". This object holds 3 fields:

pop: remove the first element (from the top).
top: access the first element (from the top).
push: add an element to the top of the Stack.
}

\author{
R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}


\seealso{
\code{\link{colShuffle}, \link{colVars}, \link{colmeans}
}
}

\examples{
x<-Stack(10,type=integer)

x$push(5)
x$push(10)
x$top()  == 10
x$pop()
x$top() == 5

y<-rnorm(10)
x<-Stack(x)

x$push(5) # length increased to 11
x$top() # access the last element that pushed, 5
x$pop() # pop the last element that pushed
}

\keyword{ Stack }
\keyword{ Environment }
