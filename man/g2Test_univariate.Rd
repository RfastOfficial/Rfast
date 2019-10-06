\name{Matrix with G-square tests of indepedence}
\alias{g2Test_univariate}
\alias{g2Test_univariate_perm}
\alias{chi2Test_univariate}

\title{
Matrix with G-square tests of indepdence
}

\description{
Matrix with G-square tests of indepdence with and without permutations.
}

\usage{
g2Test_univariate(data, dc)
g2Test_univariate_perm(data, dc, nperm) 
chi2Test_univariate(data, dc)
}

\arguments{
\item{data}{
A numerical matrix with the data. \bold{The minimum must be 0, otherwise the function can crash or will produce 
wrong results}. The data must be consecutive numbers.
}
\item{dc}{
A numerical value equal to the number of variables (or columns of the data matrix) indicating the number of distinct, unique values (or levels) of each variable.
Make sure you give the correct numbers here, otherwise the degrees of freedom will be wrong.
}
\item{nperm}{
The number of permutations. The permutations test is slower than without permutations and should be used with small sample sizes or when 
the contigency tables have zeros. When there are few variables, R's "chisq.test" function is faster, but as the number of variables increase the time difference 
with R's procedure becomes larger and larger.
}
}

\details{
The function does all the pairwise \eqn{G^2} test of independence and gives the position inside the matrix. 
The user must build the associations matrix now, similarly to the correlation matrix. See the examples of how to do that. 
The p-value is not returned, we live this to the user. See the examples of how to obtain it. 
}

\value{
A list including:
\item{statistic}{
The \eqn{G^2} or \eqn{chi^2} test statistic for each pair of variables.
}
\item{pvalue}{
This is returned when you have selected the permutation based \eqn{G^2} test.
}
\item{x}{
The row or variable of the data.
}
\item{y}{
The column or variable of the data.
}
\item{df}{
The degrees of freedom of each test.
}
}

\references{
Tsagris M. (2017). Conditional independence test for categorical data using Poisson log-linear model. 
Journal of Data Science, 15(2):347-356.

Tsamardinos, I., & Borboudakis, G. (2010). Permutation testing improves Bayesian network learning. 
In Joint European Conference on Machine Learning and Knowledge Discovery in Databases (pp. 322-337). Springer Berlin Heidelberg
}

\author{
Giorgos Borboudakis. The permutation version used a C++ code by John Burkardt.

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{g2Test}, \link{g2Test_perm}, \link{correls}, \link{univglms}
}
}

\examples{
nvalues <- 3
nvars <- 10
nsamples <- 2000
data <- matrix( sample( 0:(nvalues - 1), nvars * nsamples, replace = TRUE ), nsamples, nvars )
dc <- rep(nvalues, nvars)
system.time( g2Test_univariate(data = data, dc = dc) )
a <- g2Test_univariate(data = data, dc = dc)
pval <- pchisq(a$statistic, a$df, lower.tail = FALSE)

g <- matrix(0, nvars, nvars)
g[ cbind(a$x, a$y) ] <- a$statistic
g <- g + t(g)
diag(g) <- 0
## g ## matrix of G^2 test statistics

g<-a<-dc<-data<-NULL

}

\keyword{ G^2 test of independence }



