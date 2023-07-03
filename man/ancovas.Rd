\name{Many ANCOVAs}
\alias{ancovas}
\title{
Many ANCOVAs
}

\description{
Many ANCOVAs.
}

\usage{
ancovas(y, ina, x, logged = FALSE)
}

\arguments{
\item{y}{
A matrix with the data, where the rows denote the observations and the columns are the variables.
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
Many Analysis of covariance tests are performed. No interaction between the factor and the covariate 
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
\code{\link{ftests}, \link{ttests}, \link{anovas}
}
}

\examples{
## 100 variables, hence 100 F-tests will be performed
y <- matrix( rnorm(90 * 100), ncol = 100 )
ina <- rbinom(90, 2, 0.5) + 1
x <- rnorm(90)
a <- ancovas(y, ina, x)
\donttest{
m1 <- lm(y[, 15] ~ factor(ina) + x)
m2 <- lm(y[, 15] ~ x + factor(ina))
res<-anova(m1)
res<-anova(m2)
y <- NULL
a[15, ]  ## the same with the m2 model, but not the m1
}
}

\keyword{ Analysis of covariance }
