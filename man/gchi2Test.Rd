\name{Chi-square and G-square tests of (unconditional) indepdence}
\alias{gchi2Test}

\title{
Chi-square and G-square tests of (unconditional) indepdence
}

\description{Chi-square and G-square tests of (unconditional) indepdence.
}

\usage{
gchi2Test(x, y, logged = FALSE)

}

\arguments{
\item{x}{
A numerical vector or a factor variable with data. The data must be consecutive numbers.
}
\item{y}{
A numerical vector or a factor variable with data. The data must be consecutive numbers.
}
\item{logged}{
Should the p-values be returned (FALSE) or their logarithm (TRUE)? 
}
}

\details{
The function calculates the test statistic of the \eqn{\chi^2} and the \eqn{G^2} tests of unconditional 
independence between x and y. x and y need not be numerical vectors like in \code{\link{g2Test}}. This 
function is more close to the spirit of MASS' \code{\link{loglm}} function which calculates both statistics 
using Poisson log-linear models (Tsagris, 2017). 
}

\value{
A matrix with two rows. In each row the X2 or G2 test statistic, its p-value and the degrees of freedom are returned.
}

\references{
Tsagris M. (2017). Conditional independence test for categorical data using Poisson log-linear model. 
Journal of Data Science, 15(2):347-356. 
}

\author{
Manos Papadakis and Michail Tsagris

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com> 
and Michail Tsagris <mtsagris@uoc.gr>.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{g2Test_univariate}, \link{g2Test_univariate_perm}, \link{g2Test}
}
}

\examples{
nvalues <- 3
nvars <- 2
nsamples <- 5000
data <- matrix( sample( 0:(nvalues - 1), nvars * nsamples, replace = TRUE ), nsamples, nvars )

res<-gchi2Test(data[, 1], data[, 2])
res<-g2Test_univariate( data, rep(3, 2) )  ## G^2 test
res<-chisq.test(data[, 1], data[, 2])  ## X^2 test from R
  
data<-NULL
}

\keyword{ G^2 test of conditional independence }

