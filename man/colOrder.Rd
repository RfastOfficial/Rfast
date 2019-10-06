\name{Column and row-wise Order - Sort Indices}
\alias{colOrder}
\alias{rowOrder}
\alias{Order}
\title{
Column and row-wise Order - Sort Indices
}

\description{
Column and row-wise Order - Sort Indices.
}

\usage{
colOrder(x,stable=FALSE,descending=FALSE, parallel = FALSE)
rowOrder(x,stable=FALSE,descending=FALSE, parallel = FALSE)
Order(x,stable=FALSE,descending=FALSE,partial = NULL)
}

\arguments{
\item{x}{
A matrix with numbers or a numeric/character vector.
}
\item{stable}{
A boolean value for using a stable sorting algorithm. 
}
\item{descending}{
A boolean value (TRUE/FALSE) for sorting the vector in descending order. By default sorts the vector in ascending.
}
\item{parallel}{
A boolean value for parallel version.
}
\item{partial}{
A boolean value for partial sorting.
}
}

\details{
The function applies "order" in a column or row-wise fashion or Order a vector. If you want the same results as R's, then set "stable=TRUE" because "stable=FALSE" uses a sorting algorithm that it is not stable like R's sort. But it is faster to use the default. This verion is faster for large data, more than 300.
}

\value{
For "colOrder" and "rowOrder" a matrix with integer numbers. The result is the same as apply(x, 2, order) or apply(x, 1, order).

For "Order" sort the vector and returns the indices of each element that it has before the sorting.
The result is the same as order(x) but for the same exactly results set argument "stable" to "TRUE".
}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}

\seealso{
\code{ \link{colsums}, \link{coldiffs}, \link{colMedians}, \link{colprods}
}
}

\examples{
x <- matrix( runif(10 * 10), ncol = 10 )
res<-colOrder(x)
res<-apply(x, 2, order)
res<-rowOrder(x)
t(apply(x, 1, order))

y <- rnorm(100)
b <- Order(y)
a <- order(y)
all.equal(a,b) ## false because it is not stable
b <- Order(y,stable=TRUE)
all.equal(a,b) ## true because it is stable

x<-y<-b<-a<-NULL
}

\keyword{ Orderings }
\keyword{ matrix }
