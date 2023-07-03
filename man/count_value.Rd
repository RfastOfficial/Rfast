\name{Row - Wise matrix/vector count the frequency of a value }
\alias{colCountValues}
\alias{rowCountValues}
\alias{count_value}
\title{
Row - Wise matrix/vector count the frequency of a value
}

\description{
Row - Wise matrix/vector count the frequency of a value. 
}

\usage{
count_value(x, value)
colCountValues(x, values, parallel = FALSE)
rowCountValues(x, values, parallel = FALSE)
}

\arguments{
\item{x}{
A vector with the data (numeric or character) or a numeric matrix.
}
\item{value}{
The value, numeric or character, to check its frequency in the vector "x". 
}
\item{values}{
a vector with the values to check its frequency in the matrix "x" by row or column.
}
\item{parallel}{
Do you want to do it in parallel in C++? TRUE or FALSE. Works with every other argument.
}
}

\details{
The functions is written in C++ in order to be as fast as possible. The "x" and "value" must have the same type. The type can be numeric or character.
}

\value{
The frequency of a value/values in a vector in linear time or by row/column in a matrix.
}


\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}


\seealso{
\code{\link{Median}, \link{binary_search}, \link{Order}, \link{nth} 
}
}

\examples{
x <- rnorm(100)
value <- x[50]
count_value(x,value)
y <- sample(letters,replace=TRUE)
value <- "r"
count_value(y,value)
values <- sample(x,100,replace=TRUE)
x <- matrix(x,100,100)
res<-colCountValues(x,values)
res<-rowCountValues(x,values)
x<-value<-values<-y<-NULL
}

\keyword{ Row - Wise matrix/vector count the frequency of a value }