\name{Many one sample goodness of fit tests for categorical data}
\alias{cat.goftests}
\title{
Many one sample goodness of fit tests for categorical data
}

\description{
Many one sample goodness of fit tests for categorical data.
}

\usage{
cat.goftests(x, props, type = "gsquare", logged = FALSE)
}

\arguments{
\item{x}{
A matrix with the data, where the rows denote the samples and the columns are the variables. The data must be 
integers and be of the form 1, 2, 3, and so on. The minimum must be 1, and not zero. 
}
\item{props}{
The assumed distribution of the data. A vector or percentages summing to 1. 
}
\item{type}{
Either Pearson's \eqn{\chi^2} test ("chisquare") is used or the \eqn{G^2} test ("qsquare", default value). 
}
\item{logged}{
Should the p-values be returned (FALSE) or their logarithm (TRUE)?
}
}

\details{
Given a matrix of integers, where each column refers to a sample, the values of a categorical variable the 
function tests wether these values can be assumed to fit a specific distribution. 
}

\value{
A matrix with the test statistic and the p-value of each test.
}

%\references{
%}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{ttests}, \link{ttest}, \link{ftests}
}
}

\examples{
x <- matrix( rbinom(300 * 100, 4, 0.6), ncol = 100 ) + 1
props <- dbinom(0:4, 4, 0.6)
## can we assume that each column comes from a distribution  whose mass is given by props?
system.time( cat.goftests(x, props) )
a1 <- cat.goftests(x, props)  ## G-square test
a2 <- cat.goftests(x, props, type = "chisq")  ## Chi-square test
cor(a1, a2)
mean( abs(a1 - a2) ) 
x <- NULL
}

