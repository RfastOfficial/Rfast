\name{Hypothesis test for the distance correlation}
\alias{dcor.ttest}
\title{
Hypothesis test for the distance correlation
}

\description{
Hypothesis test for the distance correlation.
}

\usage{
dcor.ttest(x, y, logged = FALSE)
}

\arguments{
\item{x}{
A numerical matrix.
}
\item{y}{
A numerical matrix.
}
\item{logged}{
Do you want the logarithm of the p-value to be returned? If yes, set this to TRUE.
}
}

\details{
The bias corrected distance correlation is used. The hypothesis test is whether the two matrices are independent or
not. Note, that this test is size correct as both the sample size and the dimensionality goes to infinity. It will
not have the correct type I error for univariate data or for matrices with just a couple of variables.
}

\value{
A vector with 4 elements, the bias corrected distance correlation, the degrees of freedom, the test statistic and
its associated p-value. 
}

\references{
G.J. Szekely, M.L. Rizzo and N. K. Bakirov (2007). Measuring and Testing Independence 
by Correlation of Distances. Annals of Statistics, 35(6):2769-2794.
}

\author{
Manos Papadakis

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{bcdcor}, \link{dcov}, \link{edist}
}
}

\examples{
x <- as.matrix(iris[1:50, 1:4])
y <- as.matrix(iris[51:100, 1:4])
res<-dcor.ttest(x, y)
}

