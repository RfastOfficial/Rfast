\name{G-square and Chi-square test of conditional indepdence}
\alias{g2Test}
\alias{g2Test_perm}
\alias{chi2Test}

\title{
G-square test of conditional indepdence
}

\description{G-square test of conditional indepdence with and without permutations.
}

\usage{
g2tests(data, x, y, cs, dc)
g2test2_perm(data, x, y, cs, dc, nperm)
chi2test(data, x, y, cs, dc)
}

\arguments{
\item{data}{
A numerical matrix with the data. \bold{The minimum must be 0, otherwise the function can crash or will produce 
wrong results}. The data must be consecutive numbers.
}
\item{x}{
A number between 1 and the number of columns of data. This indicates which variable to take. 
}
\item{y}{
A number between 1 and the number of columns of data (other than x). This indicates the other variable whose independence with x is to be tested. 
}
\item{cs}{
A vector with the indices of the variables to condition upon. It must be non zero and between 1 and the number of variables. 
If you want unconditional independence test see \code{\link{g2Test_univariate}} and \code{\link{g2Test_univariate_perm}}.
If there is an overlap between x, y and cs you will get 0 as the value of the test statistic.
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
The functions calculates the test statistic of the \eqn{G^2} test of conditional independence between x and y conditional on a set of variable(s) cs.
}

\value{
A list including:
\item{statistic}{
The \eqn{G^2} or \eqn{chi^2} test statistic.
}
\item{df}{
The degrees of freedom of the test statistic.
}
\item{x}{
The row or variable of the data.
}
\item{y}{
The column or variable of the data.
}
}

\references{
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
\code{ \link{g2Test_univariate}, \link{g2Test_univariate_perm}, \link{correls}, \link{univglms}
}
}

\examples{
nvalues <- 3
nvars <- 10
nsamples <- 5000
data <- matrix( sample( 0:(nvalues - 1), nvars * nsamples, replace = TRUE ), nsamples, nvars )
dc <- rep(nvalues, nvars)

res<-g2Test( data, 1, 2, 3, c(3, 3, 3) )
res<-g2Test_perm( data, 1, 2, 3, c(3, 3, 3), 1000 )

dc<-data<-NULL
}

\keyword{ G^2 test of conditional independence }

