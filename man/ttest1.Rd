\name{One sample t-test for a vector}
\alias{ttest1}
\title{
One sample t-test for a vector
}

\description{
One sample t-test for a vector.
}

\usage{
ttest1(x, m, alternative = "unequal", logged = FALSE, conf = NULL)
}

\arguments{
\item{x}{
A numerical vector with the data. 
}
\item{m}{
The mean value under the null hypothesis.
}
\item{alternative}{
The alternative hypothesis, "unequal", "greater" or "less".
}
\item{logged}{
Should the p-values be returned (FALSE) or their logarithm (TRUE)?
}
\item{conf}{
If you want a confidence interval supply the confidence level.
}
}

\details{
The usual one sample t-test is implemented, only faster.
}

\value{
A list including:
\item{res}{
A two valued vector with the test statistic and its (logged) p-value.
}
\item{ci}{
In the case you supplied a number in the input argument "conf" the relevant confidence interval will 
be returned as well.
}
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> 
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{ttest}, \link{anova1}, \link{ttests}
}
}

\examples{
x = rnorm(500)
res<-t.test(x, mu = 0)
res<-ttest1(x, 0, conf = 0.95)
}

\keyword{ One sample t-test }

