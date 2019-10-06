\name{Many non parametric multi-sample tests}
\alias{kruskaltests}
\alias{cqtests}
\title{
Many multi-sample tests
}

\description{
Many multi-sample tests.
}

\usage{
kruskaltests(x, ina, logged = FALSE) 
cqtests(x, treat, block, logged = FALSE)
}

\arguments{
\item{x}{
A matrix with the data, where the rows denote the samples (and the two groups) and the columns are the variables.
}
\item{ina}{
A numerical vector with 1s, 2s, 3s and so one indicating the two groups. Be careful, the function is desinged to 
accept numbers greater than zero.
}
\item{treat}{
In the case of the Cochran's Q test, this argument plays the role of the "ina" argument.
}
\item{block}{This item denotes the subjects which are the same. Similarly to "ina" a numeric vector with 1s, 
2s, 3s and so on.
}
\item{logged}{
Should the p-values be returned (FALSE) or their logarithm (TRUE)?
}
}

\details{

The "kruskaltests" performs the Kruskal-Wallis non parametric alternative to analysis of variance test.
The "cqtests" performs the Cocrhan's Q test for the equality of more than two groups whose values are 
strictly binary (0 or 1). This is a generalisation of the McNemar's test in the multi-sample case.
}

\value{
A matrix with the test statistic and the p-value of each test.
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
\code{\link{block.anovas}, \link{ftests}
}
}

\examples{
x <- matrix( rexp(300 * 200), ncol = 200 )
ina <- rbinom(300, 3, 0.6) + 1   
system.time( kruskaltests(x, ina) ) 
x <- matrix( rbinom(300 * 200, 1, 0.6), ncol = 200 )
treat <- rep(1:3, each = 100)
block <- rep(1:3, 100)  
system.time( cqtests(x, treat, block) ) 
x <- NULL
}

\keyword{ analysis of variance }
\keyword{ Cochran's Q test }
\keyword{ non parametric statistics }