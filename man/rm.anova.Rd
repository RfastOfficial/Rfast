\name{Repeated measures anova}
\alias{rm.anova}
\title{
Repeated measures anova
}

\description{
Repeated measures anova.
}

\usage{
rm.anova(y, logged = FALSE)
}

\arguments{
\item{y}{
A matrix with the data, where each column refers to a different measurement. The rows denote the subjects.
}
\item{logged}{
Should the p-values be returned (FALSE) or their logarithm (TRUE)?
}
}

\details{
Found in Davis (2002) is the usual repeated measures ANOVA. In this case,
suppose you have taken measurements on one or more variables from the same group of people. See the example below
on how to put such data.
}

\value{
A vector with the test statistic (t-test) and its associated p-value.
}

\references{
Charles S. Davis (2002). Statistical methods for the analysis of repeated measures. Springer-Verlag, New York.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{rm.anovas}, \link{rint.reg}, \link{varcomps.mle}
}
}

\examples{
y <- c(74.5,81.5,83.6,68.6,73.1,79.4,
75.5,84.6,70.6,87.3,73.0,75.0,
68.9,71.6,55.9,61.9,60.5,61.8,
57.0,61.3,54.1,59.2,56.6,58.8,
78.3,84.9,64.0,62.2,60.1,78.7,
54.0,62.8,63.0,58.0,56.0,51.5,
72.5,68.3,67.8,71.5,65.0,67.7,
80.8,89.9,83.2,83.0,85.7,79.6)
y <- matrix(y, ncol = 6, byrow = TRUE)
res<-rm.anova(y)
}

\keyword{ Repeated measures }
\keyword{ univariate approach }
\keyword{ summary statistics }
\keyword{ regression }
