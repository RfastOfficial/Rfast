\name{Matrix with all pairs of t-tests}
\alias{allttests}
\alias{ttests.pairs}

\title{
Matrix with all pairs of t-tests
}

\description{Matrix with all pairs of t-tests.
}

\usage{
allttests(x, y = NULL, ina, logged = FALSE) 
ttests.pairs(x, logged = FALSE) 
}

\arguments{
\item{x}{
A numerical matrix with the data.
}
\item{y}{
For the case of "all.tests", if you have the second group or sample provide it here, otherwise leave it NULL.
For the case of "ttests.pairs" this is not required.
}
\item{ina}{
If you have the data in one matric then provide this indicator variable separating the samples. 
This numerical vector must contain 1s and 2s only as values. For the case of "ttests.pairs" this is not required.
}
\item{logged}{
Should the p-values be returned (FALSE) or their logarithm (TRUE)?
}
}

\details{
The function does all the pairwise t-tests assuming unequal variances (Welch's t-test). The "all.ttests" does all the pairs
formed by "cutting" the matrices x and y in two and everything between them. The "ttests.pairs" accepts a matrix x
and does all the pairs of t-tests. This is similar to the correlation matrix style. 
}

\value{
A list including:
\item{stat}{
A matrix with t-test statistic for each pair of variables.
}
\item{pvalue}{
A matrix with the corresponding p-values.
}
\item{dof}{
A matrix with the relevant degrees of freedom.
}
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
\code{\link{ttests}, \link{ftests}, \link{ttest}, \link{g2Test_univariate}
}
}

\examples{
x <- as.matrix( iris[1:100, 1:4] )
ina <- as.numeric(iris[1:100, 5])
a <- allttests(x, ina = ina)
b <- ttests.pairs(x)  ## less tests
}

\keyword{ t-tests }



