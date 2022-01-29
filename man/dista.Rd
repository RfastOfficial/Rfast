\name{Distance between vectors and a matrix}
\alias{dista}
\title{Distance between vectors and a matrix
}

\description{
Distance between vectors and a matrix.
}

\usage{
dista(xnew, x, type = "euclidean", k = 0, index = FALSE, trans = TRUE, square = FALSE)
}

\arguments{
\item{xnew}{
A matrix with some data or a vector.
}
\item{x}{
A matrix with the data, where rows denotes observations (vectors) and the columns contain the variables.
}
\item{type}{
This can be either "euclidean" or "manhattan".
}
\item{k}{
Should the k smaller distances or their indices be returned? If k > 0 this will happen.
}
\item{index}{
In case k is greater than 0, you have the option to get the indices of the k smallest distances.
}
\item{trans}{
Do you want the returned matrix to be transposed? TRUE or FALSE.
}
\item{square}{
If you choose "euclidean" as the method, then you can have the optino to return the 
squared Euclidean distances by setting this argument to TRUE. 
}
}

\details{
The target of this function is to calculate the distances between xnew and x without having to calculate the whole 
distance matrix of xnew and x. The latter does extra calculations, which can be avoided.
}

\value{
A matrix with the distances of each xnew from each vector of x. The number of rows of the xnew and and the number 
of columns of xnew are the dimensions of this matrix.
}

%\references{
%Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for compositional data. 
%In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain.
%}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm95@gmail.com>.
}


\seealso{
\code{\link{mahala}, \link{Dist}, \link{total.dist}, \link{total.dista}
}
}

\examples{
xnew <- as.matrix( iris[1:10, 1:4] )
x <- as.matrix( iris[-c(1:10), 1:4] )
a <- dista(xnew, x)
b <- as.matrix( dist( rbind(xnew, x) ) )
b <- b[ 1:10, -c(1:10) ]
sum( abs(a - b) )

## see the time
x <- matrix( rnorm(1000 * 4), ncol = 4 )
system.time( dista(xnew, x) )
system.time( as.matrix( dist( rbind(xnew, x) ) ) )

x<-b<-a<-xnew<-NULL
}

\keyword{ Distances }

