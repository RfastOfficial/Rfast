\name{Distance matrix - Sum of all pairwise distances in a distance matrix}
\alias{Dist}
\alias{total.dist}
\alias{vecdist}
\title{
Distance matrix - Sum of all pairwise distances in a distance matrix
}

\description{
Distance matrix - Sum of all pairwise distances in a distance matrix.
}

\usage{
Dist(x, method = "euclidean", square = FALSE, p = 0,vector = FALSE)
total.dist(x, method = "euclidean", square = FALSE, p = 0)
vecdist(x)
}

\arguments{
\item{x}{
A matrix with data. The distances will be calculated between pairs of rows. 
In the case of \bold{vecdist} this is a vector. 
For the \bold{haversine} distance it must be a matrix with two columns, 
the first column is the latitude and the second the longitude (in radians). 
}
\item{method}{
See details for the available methods.
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
The distance matrix is compute with an extra argument for the Euclidean distances. 
The "kullback_leibler" refers to the symmetric Kullback-Leibler divergence. 

\itemize{
\item euclidean : \eqn{  \sqrt( \sum | P_i - Q_i |^2)}
\item manhattan : \eqn{  \sum | P_i - Q_i |}
\item minimum : \eqn{  \min | P_i - Q_i |}
\item maximum : \eqn{  \max | P_i - Q_i |}
\item minkowski : \eqn{  ( \sum | P_i - Q_i |^p)^(1/p)}
\item bhattacharyya : \eqn{  - ln \sum \sqrt(P_i * Q_i)}
\item hellinger : \eqn{  2 * \sqrt( 1 - \sum \sqrt(P_i * Q_i))}
\item kullback_leibler : \eqn{  \sum P_i * log(\frac{P_i}{Q_i}) + \sum P_i * log(\frac{Q_i}{P_i})}
\item jensen_shannon : \eqn{  0.5 * ( \sum P_i * log(\frac{2 * P_i}{P_i + Q_i}) + \sum Q_i * log(\frac{2 * Q_i}{P_i + Q_i})) }
\item haversine : \eqn{ 2 * R * \arcsin(\sqrt(\sin((lat_2 - lat_1)/2)^2 + \cos(lat_1) * \cos(lat_2) * \sin((lon_2 - lon_1)/2)^2)) }
\item canberra : \eqn{  \sum | P_i - Q_i | / (P_i + Q_i)}
\item chi_square \eqn{X}^2 : \eqn{  \sum\frac{(P_i - Q_i )^2}{P_i + Q_i}}
\item soergel : \eqn{  \sum | P_i - Q_i | / \sum \max(P_i , Q_i)}
\item sorensen : \eqn{  \sum | P_i - Q_i | / \sum (P_i + Q_i)}
\item cosine : \eqn{ \sum (P_i * Q_i) / \sqrt(\sum P_i^2) * \sqrt(\sum Q_i^2)}
\item wave_hedges : \eqn{  \sum | P_i - Q_i | / \max(P_i , Q_i)}
\item motyka : \eqn{  \sum \min(P_i , Q_i) / (P_i + Q_i)}
\item harmonic_mean : \eqn{ 2 * \sum (P_i * Q_i) / (P_i + Q_i)}
\item jeffries_matusita : \eqn{  \sqrt( 2 - 2 * \sum \sqrt(P_i * Q_i))}
\item gower : \eqn{  1/d * \sum | P_i - Q_i |}
\item kulczynski : \eqn{  1 / \sum | P_i - Q_i | / \sum \min(P_i , Q_i)}
}
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

\keyword{ Distance matrix - Sum of all pairwise distances in a distance matrix }
\keyword{ Euclidean distance }
\keyword{ Manhattan distance }
\keyword{ Canberra distance }
\keyword{ Hellinger distance }