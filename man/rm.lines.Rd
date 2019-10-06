\name{Many regression based tests for single sample repeated measures}
\alias{rm.lines}
\alias{rm.anovas}
\title{
Many regression based tests for single sample repeated measures
}

\description{
Many regression based tests for single sample repeated measures.
}

\usage{
rm.lines(y, x, logged = FALSE)
rm.anovas(y, x, logged = FALSE)
}

\arguments{
\item{y}{
A matrix with the data, where each column refers to a different sample of subjects. For example, the first
column is the repeated measurements of a sample of subjects, the second column contains repeated measurements of
a second sample of subjects and so on. Within each column, the measurements of each subjects are stacked one
upon the other. Say for examples there are n subjects and each of them has been measured d times (in time or at
different experimental conditions). We put these in a matrix with just one column. The first d rows are the
measurements of subject 1, the next d rows are the measurements of subject 2 and so on.
}
\item{x}{
A numerical vector with time (usually) or the the predictor variable. For example the temperature,
or the pressure. See the details for more information. Its length is equal to the time points for example,
i.e. it must not have the same length as the number of rows of y. For the "rm.lines" this is a continuous
variable.

For the "rm.anovas" this is treated as a categorical variable, indicating say the type of experimental condition,
but no difference between the points is important. Hence, for this function only, x can also be a facto variable.
}
\item{logged}{
Should the p-values be returned (FALSE) or their logarithm (TRUE)?
}
}

\details{
In order to see whether the repeated measurements are associated with a single covariate, e.g. time we
perform many regressions and each time calculate the slope. For each subject, its regression slope with
the covariate is calculated. In the end a t-test for the hypothesis that the average slopes is zero is performed.
The regression slopes ignore that the measurements are not independent, but note that the slopes are independent,
because they come from different subjects. This is a simple, summary statistics based approach found in Davis
(2002), yet it can provide satisfactory results.

The second approach ("rm.anovas") found in Davis (2002) is the usual repeated measures ANOVA. In this case,
suppose you have taken measurements on one or more variables from the same group of people. See the example below
on how to put such data.
}

\value{
A matrix with the test statistic (t-test) and its associated p-value.
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
\code{ \link{rint.regbx}, \link{rint.reg}, \link{varcomps.mle}
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
y <- as.matrix(y)
### the first 6 measurements are from subject 1, measurments 7-12 are from subject 2,
## measurements 13-18 are from subject 3 and so on.
x <- c(-10, 25, 37, 50, 65, 80) ## all subjects were measured at the same time points
res<-rm.lines(y, x) ## Is linear trend between the measurements and the temperature?
res<-rm.anovas(y, x)  ## Tests whether the means of the individuals are the same
## the temperature is treated as categorical variable here.

## fake example
y <- matrnorm(10, 4)
## the y matrix contains 4 repeated measurements for each of the 10 persons.
x <- 1:4
## we stack the measurements of each subject, one under the other in a matrix form.
y1 <- matrix( t(y) )
res<-rm.anovas(y1, x)  ## perform the test
z <- matrix( rnorm(20 * 8), ncol = 2) ## same example, but with 2 sets of measurements.
res<-rm.anovas(z, x)
}

\keyword{ Repeated measures }
\keyword{ univariate approach }
\keyword{ summary statistics }
\keyword{ regression }
