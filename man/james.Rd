\name{James multivariate version of the t-test}
\alias{james}
\title{
James multivariate version of the t-test
}

\description{
James test for testing the equality of two population mean vectors without assuming equality of the covariance matrices.
}

\usage{
james(y1, y2, a = 0.05, R = 1)
}

\arguments{
  \item{y1}{
  A matrix containing the Euclidean data of the first group.
}
  \item{y2}{
  A matrix containing the Euclidean data of the second group.
}
  \item{a}{
The significance level, set to 0.05 by default.
}
\item{R}{
If R is 1 the classical James test is returned. 
If R is 2 the MNV modficiation is implemented.
}
}

\details{
Multivariate analysis of variance without assuming equality of the covariance matrices. The p-value can be calculated either asymptotically or via bootstrap. The James test (1954) or a modification proposed by Krishnamoorthy and Yanping (2006) is implemented. The James test uses a corected chi-square distribution, whereas the modified version uses an F distribution.
}

\value{
A list including:
\item{note}{
A message informing the user about the test used.
}
\item{mesoi}{
The two mean vectors.
}
\item{info}{
The test statistic, the p-value, the correction factor and the corrected critical 
value of the chi-square distribution if the James test has been used or,
the test statistic, the p-value, the critical value and the degrees of freedom 
(numerator and denominator) of the F distribution if the modified James test 
has been used.
}
}

\references{
G.S. James (1954). Tests of Linear Hypothese in Univariate and Multivariate Analysis
when the Ratios of the Population Variances are Unknown. Biometrika, 41(1/2): 19-43

Krishnamoorthy K. and Yanping Xia.  On Selecting Tests for Equality of Two Normal Mean Vectors (2006).
Multivariate Behavioral Research 41(4): 533-548
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{mv.eeltest2}
}
}

\examples{
james( as.matrix(iris[1:25, 1:4]), as.matrix(iris[26:50, 1:4]), R = 1 )
james( as.matrix(iris[1:25, 1:4]), as.matrix(iris[26:50, 1:4]), R = 2 )
}

\keyword{ Multivariate analysis of variance }
\keyword{ Unequality of the covariance matrices }
