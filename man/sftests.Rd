\name{Many Shapiro-Francia normality tests}
\alias{sftests}
\alias{sftest}
\title{
Many Shapiro-Francia normality tests
}

\description{
Many Shapiro-Francia normality tests.
}

\usage{
sftests(x, logged = FALSE)
sftest(x, logged = FALSE)
}

\arguments{
\item{x}{
A matrix with the data, where the rows denote the observations and the columns are the variables. 
In the case of a single sample, then this must be a vector and "sftest" is to be used.
}
\item{logged}{
Should the p-values be returned (FALSE) or their logarithm (TRUE)?
}
}

\details{
The Shapiro-Francia univariate normality test is performed for each column (variable) of the matrix x. 
}

\value{
A matrix with the squared correlation between the ordered values and the standard normal ordered statistics, 
the test statistic and the p-value of each test. If the "sftest" has been used, the output is a vector with these three elements. 
}

\references{
Royston J. P. (1983). A simple method for evaluating the Shapiro-Francia W' test of non-normality. 
The Statistician, 32(3): 297-300.

Mbah A. K. & Paothong A. (2015). Shapiro-Francia test compared to other normality test using expected p-value. 
Journal of Statistical Computation and Simulation, 85(15): 3002-3016.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{ttests}, \link{ttest}, \link{ftests}
}
}

\examples{
x <- matrnorm(200, 100)
system.time( sftests(x) )
a <- sftests(x)
mean(a[, 3]<0.05) 
x <- rnorm(100)
res<-sftest(x)
}

\keyword{ Univariate normality test }
\keyword{ Shapiro-Francia }
