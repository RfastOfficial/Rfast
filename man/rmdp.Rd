\name{High dimensional MCD based detection of outliers}
\alias{rmdp}
\title{
High dimensional MCD based detection of outliers
}

\description{
High dimensional MCD based detection of outliers.
}

\usage{
rmdp(y, alpha = 0.05, itertime = 100)
}

\arguments{
\item{y}{
A matrix with numerical data with more columns (p) than rows (n), i.e. n<p.
}
\item{alpha}{
The significance level, i.e. used to decide whether an observation is said to be considered a possible outlier. 
The default value is 0.05.
}
\item{itertime}{
The number of iterations the algorithm will be ran. The higher the sample size, the larger this number must be. 
With 50 observations in \eqn{R^1000} maybe this has to be 1000 in order to produce stable results.  
}
}

\details{
High dimensional outliers (n<<p) are detected using a properly constructed MCD. The variances of the variables are 
used and the determinant is simply their product. 
}

\value{
A list including:
runtime = runtime, dis = dis, wei = wei
\item{runtime}{
The duration of the process.
}
\item{dis}{
The final estimated Mahalanobis type normalised distances.
}
\item{wei}{
A bollean variable vector specifying whether an observation is "clean" (TRUE) or a possible outlier (FALSE).
}
\item{cova}{
The estimated covatriance matrix.
}
}

\references{
Ro K., Zou C., Wang Z. and Yin G. (2015). Outlier detection for high-dimensional data. Biometrika, 102(3):589-599.
}

\author{
Initial R code: Changliang Zou <nk.chlzou@gmail.com>
R code modifications: Michail Tsagris <mtsagris@yahoo.gr>
C++ implementation: Manos Papadakis <papadakm95@gmail.com>
Documentation: Michail Tsagris <mtsagris@yahoo.gr> and Changliang Zhou <nk.chlzou@gmail.com>
}

\seealso{
\code{ \link{colmeans}, \link{colVars}, \link{colMedians} 
}
}

\examples{
x <- matrix(rnorm(50 * 400), ncol = 400)
a <- rmdp(x, itertime = 500)

x<-a<-NULL
}

\keyword{ High dimensional data }
\keyword{ outliers }
\keyword{ MCD estimation }