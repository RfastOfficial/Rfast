\name{Analysis of covariance}
\alias{ancova1}
\title{
Analysis of covariance
}

\description{
Analysis of covariance
}

\usage{
ancova1(y, ina, x, logged = FALSE)
}

\arguments{
\item{y}{
A numerical vector with the data, the response variable.
}
\item{ina}{
A numerical vector with 1s, 2s, 3s and so one indicating the two groups. Be careful, the function 
is desinged to accept numbers greater than zero. 
}
\item{x}{
A numerical vector whose length is equal to the number of rows of y. This is the covariate. 
}
\item{logged}{
Should the p-values be returned (FALSE) or their logarithm (TRUE)?
}
}

\details{
Analysis of covariance is performed. No interaction between the factor and the covariate 
is tested. Only the main effects. The design need not be balanced. The values of ina need not have the 
same frequency. The sums of squares have been adjusted to accept balanced and unbalanced designs. 
}

\value{
A matrix with the test statistic and the p-value for the factor variable and the covariate. 
}

\references{
D.C. Montgomery (2001). Design and analysis of experiments (5th Edition). New York: John Wiley & Sons 
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{ancovas}, \link{ftests}, \link{ttests}, \link{anova1}
}
}

\examples{
y <- rnorm(90)
ina <- rbinom(90, 2, 0.5) + 1
x <- rnorm(90)
system.time( a <- ancova1(y, ina, x) )
}

\keyword{ Analysis of covariance }
