\name{Energy distance between matrices}
\alias{edist}
\title{
Energy distance between matrices
}

\description{
Energy distance between matrices.
}

\usage{
edist(x, y = NULL)
}

\arguments{
\item{x}{
A matrix with numbers or a list with matrices. 
}
\item{y}{
A second matrix with data. The number of columns of x and y must match. The number of rows can be different. 
}
}

\details{
This calculates the energy distance between two matrices. It will work even for tens of thousands of rows, 
it will just take some time. See the references for more information. If you have many matrices and want to 
calculate the distance matrix, then put them in a list and use the function.
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

Tsagris M. and Papadakis M. (2025). Fast and light-weight energy statistics using the R package Rfast.
https://arxiv.org/abs/2501.02849
}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}

\seealso{
\code{ \link{dcov}, \link{Dist}, \link{dista}
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
