\name{Many 2 sample tests}
\alias{ttests}
\alias{mcnemars}
\alias{var2tests}
\title{
Many 2 sample tests tests
}

\description{
It performs very many 2 sample tests.
}

\usage{
ttests(x, y = NULL, ina, alternative = "unequal", paired = FALSE, logged = FALSE, parallel = FALSE)
mcnemars(x, y = NULL, ina, logged = FALSE) 
var2tests(x, y = NULL, ina, alternative = "unequal", logged = FALSE) 
}

\arguments{
\item{x}{
A matrix with the data, where the rows denote the samples and the columns are the variables.
}
\item{y}{
A second matrix with the data of the second group. If this is NULL (default value) then the argument ina must be supplied. Notice that when you supply the two matrices the procedure is two times faster.
} 
\item{ina}{
A numerical vector with 1s and 2s indicating the two groups. Be careful, the function is designed to accept only these two numbers. In addition, if your "y" is NULL, you must specify "ina".
}
\item{alternative}{
The type of hypothesis to be checked, "equal", "greater", "less". 
}
\item{paired}{
If the groups are not independent paired t-tests should be performed and this must be TRUE, otherwise, 
leave it FALSE. In this case, the two groups must have equal smaple sizes, otherwise no test will be performed.
}
\item{logged}{
Should the p-values be returned (FALSE) or their logarithm (TRUE)?
}
\item{parallel}{
Should parallel implentations take place in C++? The default value is FALSE.
}
}

\details{
For the ttests, if the groups are independent, the Welch's t-test (without assuming equal variances) is 
performed. Otherwise many paired t-tests are performed. The McNemar's test requires a number of observations, 
at least 30 would be good in order for the test to have some power and be size corect.  
}

\value{
A matrix with the test statistic, the degrees of freedom (if the groups are independent) and the p-value 
(or their logarithm) of each test.
}

\references{
B. L. Welch (1951). On the comparison of several mean values: an alternative approach. Biometrika, 38(3/4), 330-336.

McNemar Q. (1947). Note on the sampling error of the difference between correlated proportions or percentages. Psychometrika. 12(2):153-157. 
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{ftests}, \link{anovas}, \link{ttest}
}
}

\examples{
## 1000 variables, hence 1000 t-tests will be performed
x = matrnorm(100, 100)
## 100 observations in total
ina = rbinom(100, 1, 0.6) + 1   ## independent samples t-test
ttests(x, ina = ina)
x1 = x[ina == 1, ]
x2 = x[ina == 2, ]
ttests(x1, x2)
x <- NULL
}

\keyword{ t-tests }
\keyword{ McNemar's test }
\keyword{ 2 variances test }

