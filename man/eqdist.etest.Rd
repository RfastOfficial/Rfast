\name{Energy test of equal univariate distributions}
\alias{eqdist.etest}
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
A numerical vector or a numerical matrix. 
}
\item{x}{
A numerical vector or a numerical matrix. 
}
\item{R}{
The number of permutations to perform.
}
}

\details{
The test performs the energy test of equal univariate distributions and the p-value is computed via permutations.
Both the univariate and multivariate cases are memory-saving, the univariate case is pretty fast, but the 
multivariate case is not fast enough.
}

\value{
The permutation based p-value.
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
https://www.researchgate.net/publication/387583091_Fast_and_light-weight_energy_statistics_using_the_R_package_Rfast
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>.
}

\seealso{
\code{ \link{edist}, \link{dcov}
}
}

\examples{
y <- iris[1:50, 1] 
x <- iris[1:50, 2] 
eqdist.etest(y, x, R = 499)
}

\keyword{ Energy distances }
\keyword{ total sum }
