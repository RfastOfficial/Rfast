\name{Odds ratio and relative risk}
\alias{odds.ratio}
\alias{rel.risk}
\title{
Odds ratio and relative risk
}

\description{
Odds ratio and relative risk.
}

\usage{
odds.ratio(x, a = 0.05, logged = FALSE)
rel.risk(x, a = 0.05, logged = FALSE)
}

\arguments{
\item{x}{
A 2 x 2 matrix or a vector with 4 elements. In the case of the vector make sure it corresponds to the correct table.
}
\item{a}{
The significance level, set to 0.05 by default.
}
\item{logged}{
Should the p-values be returned (FALSE) or their logarithm (TRUE)?
}
}

\details{
The odds ratio and the confidence interval are calculated.
}

\value{
A list including:
\item{res}{
The estimated odds ratio and the p-value for the null hypothesis test that it is equal to 1.
}
\item{ci}{
The (1-a)\% confidence interval for the true value of the odds ratio.
}
}

\references{
Mosteller Frederick (1968). Association and Estimation in Contingency Tables. Journal of the American Statistical 
Association. 63(321):1-28.
 
Edwards A.W.F. (1963). The measure of association in a 2x2 table. Journal of the Royal Statistical Society, Series A. 
126(1):109-114. 
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{odds}, \link{g2Test}
}
}

\examples{
x <- rpois(4, 30)+2
res<-odds.ratio(x)
res<-odds.ratio( matrix(x, ncol = 2) )
}

\keyword{ Odds ratio }

