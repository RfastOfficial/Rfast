\name{Many 2 sample proportions tests}
\alias{proptests}
\title{
Many 2 sample proportions tests
}

\description{
It performs very many 2 sample proportions tests.
}

\usage{
proptests(x1, x2, n1, n2) 
}

\arguments{
\item{x1}{
A vector with the successes of the one group.
}
\item{x2}{
A vector with the successes of the one group.
}
\item{n1}{
A vector with the number of trials of the one group.
}
\item{n2}{
A vector with the number of trials of the one group.
}
}

\details{
The 2-sample proportions test is performed for each pair of proportions of teh two groups. 
}

\value{
A matrix with the proportions of each group (two columns), the test statistic and the p-value of each test.
}

\references{
B. L. Welch (1951). On the comparison of several mean values: an alternative approach. Biometrika, 38(3/4), 330-336.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{ttests}, \link{ftests}, \link{colVars}
}
}

\examples{
## 10000 variables, hence 10000 t-tests will be performed
set.seed(12345)
x1 <- rpois(500, 5)
x2 <- rpois(500, 5)
n1 <- rpois(1000, 40)
n2 <- rpois(1000, 40)
a <- proptests(x1, x2, n1, n2)
mean(a[, 4]<0.05)

x1 <- rbinom(500, 500, 0.6)
x2 <- rbinom(500, 500, 0.6)
b <- proptests(x1, x2, 500, 500)
mean(b[, 4]<0.05)
}

\keyword{ 2 sample proportions tests }

