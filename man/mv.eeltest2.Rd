\name{Exponential empirical likelihood hypothesis testing for two mean vectors}
\alias{mv.eeltest2}
\title{
Exponential empirical likelihood hypothesis testing for two mean vectors
}

\description{
Exponential empirical likelihood hypothesis testing for two mean vectors.
}

\usage{
mv.eeltest2(y1, y2, tol = 1e-07, R = 0)
}

\arguments{
\item{y1}{
A matrix containing the Euclidean data of the first group.
}
\item{y2}{
A matrix containing the Euclidean data of the second group.
}
\item{tol}{
The tolerance level used to terminate the Newton-Raphson algorithm.
}
\item{R}{
If R is 0, the classical chi-sqaure distribution is used, if R = 1, the corrected chi-square distribution (James, 1954) is used and if R = 2, the modified F distribution (Krishnamoorthy and Yanping, 2006) is used.
}
}

\details{
Exponential empirical likelihood is a non parametric hypothesis testing procedure for one sample. The generalisation to two (or more samples) is via searching for the mean vector that minimises the sum of the two test statistics.
}

\value{
A list including:
\item{test}{
The empirical likelihood test statistic value.}
\item{modif.test}{
The modified test statistic, either via the chi-square or the F distribution.
}
\item{pvalue}{
The p-value.
}
\item{iters}{
The number of iterations required by the newton-Raphson algorithm.
}
\item{mu}{
The estimated common mean vector.
}
}

\references{
Jing Bing-Yi and Andrew TA Wood (1996). Exponential empirical likelihood is not Bartlett correctable. Annals of Statistics 24(1): 365-369.

G.S. James (1954). Tests of Linear Hypothese in Univariate and Multivariate Analysis
when the Ratios of the Population Variances are Unknown. Biometrika, 41(1/2): 19-43.

Krishnamoorthy K. and Yanping Xia (2006).  On Selecting Tests for Equality of Two Normal Mean Vectors.
Multivariate Behavioral Research 41(4): 533-548.

Owen A. B. (2001). Empirical likelihood. Chapman and Hall/CRC Press.

Amaral G.J.A., Dryden I.L. and Wood A.T.A. (2007). Pivotal bootstrap methods for k-sample problems in directional statistics and shape analysis. 
Journal of the American Statistical Association 102(478): 695-707.

Preston S.P. and Wood A.T.A. (2010). Two-Sample Bootstrap Hypothesis Tests for Three-Dimensional Labelled Landmark Data. Scandinavian Journal of 
Statistics 37(4): 568-587.

Tsagris M., Preston S. and Wood A.T.A. (2017). Nonparametric hypothesis testing for equality of means on the
simplex. Journal of Statistical Computation and Simulation, 87(2): 406-422.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr>.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{james}, \link{mv.eeltest1}
}
}

\examples{
res<-mv.eeltest2( as.matrix(iris[1:25, 1:4]), as.matrix(iris[26:50, 1:4]), R = 0 )
res<-mv.eeltest2( as.matrix(iris[1:25, 1:4]), as.matrix(iris[26:50, 1:4]), R = 1 )
}

\keyword{ Multivariate hypothesis testing }
\keyword{ non parametric test }
