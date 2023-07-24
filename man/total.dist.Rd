\name{Sum of all pairwise distances in a distance matrix}
\alias{total.dist}
\alias{total.dista}
\title{
Sum of all pairwise distances in a distance matrix
}

\description{
Sum of all pairwise distances in a distance matrix.
}

\usage{
total.dista(x, y, square = FALSE)  
}

\arguments{
\item{x}{
A matrix with numbers.
}
\item{y}{
A second matrix with data. The number of comlumns of this matrix must be the same with the matrix x. The number of rows can be different. 
}
\item{method}{
This is either "euclidean", "manhattan", "canberra1", "canberra2", "minimum", "maximum", "minkowski","bhattacharyya", "hellinger", "total_variation" or "kullback_leibler/jensen_shannon".
The last two options are basically the same. 
}
\item{square}{
If you choose "euclidean" or "hellinger" as the method, then you can have the option to return the 
squared Euclidean distances by setting this argument to TRUE. 
}
\item{p}{
This is for the the Minkowski, the power of the metric. 
}
}

\details{
In order to do the total.dist one would have to calcualte the distance matrix and sum it. We do this internally in C++ without creating the matrix. 
For the total.dista it is the same thing.
}

\value{
A numerical value, the sum of the distances. 
}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}

\seealso{
\code{ \link{Dist}, \link{dista}
}
}

\examples{
x <- matrix( rnorm(50 * 10), ncol = 10 )
res<-total.dist(x) 
y <- matrix( rnorm(40 * 10), ncol = 10)
res<-total.dista(x, y)
res<-total.dista(y, x)

x<-y<-NULL
}

\keyword{ Distances}
\keyword{ total sum }
