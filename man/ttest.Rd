\name{Many one sample tests}
\alias{ttest}
\alias{proptest}
\alias{vartest}

\title{
Many one sample tests
}

\description{
Many one sample tests.
}

\usage{
proptest(x, n, p, alternative = "unequal", logged = FALSE) 
ttest(x, m, alternative = "unequal", logged = FALSE, conf = NULL)
vartest(x, sigma, alternative = "unequal", logged = FALSE, conf = NULL) 
}

\arguments{
\item{x}{
A matrix with numerical data. Each column of the matrix corresponds to a sample, or a group. In the case of the "proptest" this is a vector integers ranging from 0 up to n. It is the number of "successes". 
}
\item{n}{
This is for the "proptest" only and is a vector with integer numbers specifying the number of tries for the proptest. Its size is equal to the size of x.
}
\item{p}{
A vector with the assumed probabilities of success in the "proptest". Its size is equal to the number of colums of the matrix x.
}
\item{m}{
A vector with the assumed means. Its size is equal to the number of colums of the matrix x.
}
\item{sigma}{
A vector with assumed variances. Its size is equal to the number of colums of the matrix x.
}
\item{alternative}{
The type of hypothesis to be checked. Equal to ("unequal"), grater than("greater") or less than ("less") the assumed parameter.
}
\item{logged}{
Should the p-values be returned (FALSE) or their logarithm (TRUE)?
}
\item{conf}{
If you want confidence intervals to be returned specify the confidence level, otherwise leave it NULL.
}
}

\details{
Despite the functions having been written in R, they are very fast.  
}

\value{
For all tests except for the "sftests" a matrix with two colums, the test statistic and the p-value respectively. 
}

%\references{
%}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{ftests}, \link{ttests}
}
}

\examples{
R <- 100
## protest
x <- rbinom(R, 50, 0.6)
n <- rep(50, R)
p <- rep(0.6, R)
a1 <- proptest(x, n, p, "unequal", logged = FALSE)
res<-sum( a1[, 2] < 0.05 ) / R

## vartest
x <- matrnorm(100, 100)
a2 <- vartest(x, rep(1, R) )
res<-sum( a2[, 2] < 0.05 )

## ttest
a4 <- ttest(x, numeric(R) )
res<-sum(a4[, 2] < 0.05) / R
x <- NULL
}

\keyword{ Hypothesis testing }
\keyword{ one sample }
\keyword{ t-test }
\keyword{ variance test }
\keyword{ proportion test }
