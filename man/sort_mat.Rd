\name{Sorting of the columns-rows of a matrix}
\alias{colSort}
\alias{rowSort}
\alias{sort_mat}
\title{
Sorting of the columns-rows of a matrix 
}

\description{
Fast sorting of the columns-rows of a matrix.
}

\usage{
colSort(x, descending = FALSE, stable = FALSE,parallel=FALSE)
rowSort(x, descending = FALSE, stable = FALSE,parallel=FALSE)
sort_mat(x,by.row=FALSE,descending=FALSE,stable=FALSE,parallel=FALSE)
}

\arguments{
\item{x}{
A numerical matrix or data.frame with data.
}
\item{descending}{
If you want the sorting in descending order, set this to TRUE.
}
\item{stable}{
If you the stable version, so that the results are the same as R's (in the case of ties) set this to TRUE.
If this is TRUE, the algorithm is a bit slower.
}
\item{parallel}{
Do you want to do it in parallel in C++? TRUE or FALSE. Works with every other argument.
}
\item{by.row}{
TRUE or FALSE for applying sort in rows or column.
}
}


\value{
The matrix with its columns-rows (or rows) independently sorted.
}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}

\seealso{
\code{\link{nth}, \link{colMaxs}, \link{colMins}, \link{colrange}, \link{sort_cor_vectors}, 
\link{sort_unique}
}
}

\examples{
x <- matrix( rnorm(100 * 200), ncol = 200 )
s1 <- colSort(x)
#s2 <- apply(x, 2, sort)
#all.equal(as.vector(s1), as.vector(s2))

x<-NULL
}

\keyword{ Sorting }


