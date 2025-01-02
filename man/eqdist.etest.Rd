\name{Energy test of equal univariate distributions}
\alias{edist}
\title{
Energy test of equal univariate distributions
}

\description{
Energy test of equal univariate distributions.
}

\usage{
eqdist.etest(y, x, R = 999)
}

\arguments{
\item{y}{
A numerical vector. 
}
\item{x}{
A numerical vector. 
}
\item{R}{
The number of permutations to perform.
}
}

\details{
The test performs the energy test of equal univariate distributions and the p-value is computed via pmerutations.
}

\value{
If "x" is matrix, a numerical value, the energy distance. If "x" is list, a matrix with all pairwsie distances of the matrices.
}

\references{
Szekely G. J. and Rizzo M. L. (2004) Testing for Equal Distributions in High Dimension, InterStat, November (5).

Szekely G. J. (2000) Technical Report 03-05, E-statistics: Energy of Statistical Samples, Department of 
Mathematics and Statistics, Bowling Green State University. 

Sejdinovic D., Sriperumbudur B., Gretton A. and Fukumizu, K. (2013). Equivalence of distance-based and
RKHS-based statistics in hypothesis testing. The Annals of Statistics, 41(5): 2263--2291.

Szekely G. J. and Rizzo M. L. (2023). The Energy of Data and Distance Correlation. 
Chapman and Hall/CRC.
}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}

\seealso{
\code{ \link{dvar}, \link{Dist}, \link{dista}
}
}

\examples{
x <- as.matrix( iris[1:50, 1:4] )
y <- as.matrix( iris[51:100, 1:4] )
res<-edist(x, y)
z <- as.matrix(iris[101:150, 1:4])
a <- list()
a[[ 1 ]] <- x
a[[ 2 ]] <- y
a[[ 3 ]] <- z
res<-edist(a)

x<-y<-z<-a<-NULL
}

\keyword{ Energy distances }
\keyword{ total sum }
