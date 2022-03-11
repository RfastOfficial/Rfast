\name{Distance matrix}
\alias{Dist}
\alias{vecdist}
\title{
Distance matrix
}

\description{
Distance matrix.
}

\usage{
Dist(x, method = "euclidean", square = FALSE, p = 0,vector = FALSE)
vecdist(x)
}

\arguments{
\item{x}{
A matrix with data. The distances will be calculated between pairs of rows. 
In the case of \bold{vecdist} this is a vector. For the haversine distance it must be 
a matrix with two columns, the first column is the latitude and the second the longitude. 
}
\item{method}{
This is either "euclidean", "manhattan", "canberra1", "canberra2", "minimum", "maximum", 
"minkowski", "bhattacharyya", "hellinger", "kullback_leibler", "jensen_shannon" or "haversine". 
}
\item{square}{
If you choose "euclidean" or "hellinger" as the method, then you can have the option to return the 
squared Euclidean distances by setting this argument to TRUE. 
}
\item{p}{
This is for the the Minkowski, the power of the metric. 
}
\item{vector}{
For return a vector instead a matrix.
}
}

\details{
The distance matrix is computer with an extra argument for the Euclidean distances. 
The "kullback_leibler" refers to the symmetric Kullback-Leibler divergence. 
}

\value{
A square matrix with the pairwise distances. 
}

\references{
Mardia K. V., Kent J. T. and Bibby J. M. (1979). Multivariate Analysis. Academic Press. 
}

\author{
Manos Papadakis.

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}


\seealso{
\code{\link{dista}, \link{colMedians}
}
}

\examples{
x <- matrix(rnorm(50 * 10), ncol = 10)
a1 <- Dist(x)
a2 <- as.matrix( dist(x) )

x<-a1<-a2<-NULL
}

\keyword{ Distance matrix }
\keyword{ Euclidean distance }
\keyword{ Manhattan distance }
\keyword{ Canberra distance }
\keyword{ Hellinger distance }