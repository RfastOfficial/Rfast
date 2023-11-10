\name{Distance between vectors and a matrix - Sum of all pairwise distances in a distance matrix.}
\alias{dista}
\alias{total.dista}
\title{Distance between vectors and a matrix - Sum of all pairwise distances in a distance matrix.
}

\description{
Distance between vectors and a matrix - Sum of all pairwise distances in a distance matrix..
}

\usage{
dista(xnew, x, type = "euclidean", k = 0, index = FALSE, 
 trans = TRUE, square = FALSE, p = 0, parallel = FALSE)
total.dista(xnew, x, type = "euclidean", k = 0,
 square = FALSE, p = 0, parallel = FALSE)  
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
If you choose "euclidean" or "hellinger" as the method, then you can have the option to return the 
squared Euclidean distances by setting this argument to TRUE. 
}
\item{p}{
This is for the the Minkowski, the power of the metric. 
}
\item{parallel}{
For methods \bold{kullback_leibler}, \bold{jensen_shannon} and \bold{itakura_saito}, you can run the 
algorithm in parallel. 
}
}

\details{
The target of this function is to calculate the distances between xnew and x without having to calculate the whole 
distance matrix of xnew and x. The latter does extra calculations, which can be avoided.

\itemize{
\item euclidean : \eqn{ \sum \sqrt( \sum | P_i - Q_i |^2)}
\item manhattan : \eqn{ \sum \sum | P_i - Q_i |}
\item minimum : \eqn{ \sum \min | P_i - Q_i |}
\item maximum : \eqn{ \sum \max | P_i - Q_i |}
\item minkowski : \eqn{ \sum ( \sum | P_i - Q_i |^p)^(1/p)}
\item bhattacharyya : \eqn{ \sum - ln \sum \sqrt(P_i * Q_i)}
\item hellinger : \eqn{ \sum 2 * \sqrt( 1 - \sum \sqrt(P_i * Q_i))}
\item kullback_leibler : \eqn{ \sum \sum P_i * log(P_i / Q_i)}
\item jensen_shannon : \eqn{ \sum 0.5 * ( \sum P_i * log(2 * P_i / P_i + Q_i) + \sum Q_i * log(2 * Q_i / P_i + Q_i))}
\item canberra : \eqn{ \sum \sum | P_i - Q_i | / (P_i + Q_i)}
\item chi_square \eqn{X}^2 : \eqn{ \sum \sum ( (P_i - Q_i )^2 / (P_i + Q_i) )}
\item soergel : \eqn{ \sum \sum | P_i - Q_i | / \sum \max(P_i , Q_i)}
\item sorensen : \eqn{ \sum \sum | P_i - Q_i | / \sum (P_i + Q_i)}
\item cosine : \eqn{ \sum (P_i * Q_i) / \sqrt(\sum P_i^2) * \sqrt(\sum Q_i^2)}
\item wave_hedges : \eqn{ \sum \sum | P_i - Q_i | / \max(P_i , Q_i)}
\item motyka : \eqn{ \sum \sum \min(P_i , Q_i) / (P_i + Q_i)}
\item harmonic_mean : \eqn{ 2 * \sum (P_i * Q_i) / (P_i + Q_i)}
\item jeffries_matusita : \eqn{ \sum \sqrt( 2 - 2 * \sum \sqrt(P_i * Q_i))}
\item gower : \eqn{ \sum 1/d * \sum | P_i - Q_i |}
\item kulczynski : \eqn{ \sum 1 / \sum | P_i - Q_i | / \sum \min(P_i , Q_i)}
}
}

\value{
A matrix with the distances of each xnew from each vector of x. The number of rows of the xnew and and the number 
of columns of xnew are the dimensions of this matrix.
}

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
dista(xnew, x)
as.matrix( dist( rbind(xnew, x) ) )

x<-b<-a<-xnew<-NULL
}

\keyword{ Distances }