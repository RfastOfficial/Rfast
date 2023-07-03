\name{Column-wise kurtosis and skewness coefficients}
\alias{colkurtosis}
\alias{colskewness}
\title{
Column-wise kurtosis and skewness coefficients
}

\description{
Column-wise kurtosis and skewness coefficients.
}

\usage{
colkurtosis(x, pvalue = FALSE)

colskewness(x, pvalue = FALSE)
}

\arguments{
\item{x}{
A matrix with the data, where the rows denote the samples and the columns are the variables.
}
\item{pvalue}{
If you want a hypothesis test that the skewness or kurtosis are significant set this to TRUE. This checks
whether the skewness is significantly different from 0 and whether the kurtosis is significantly different from
3.  
}
}

\details{
The skewness and kurtosis coefficients are calculated. For the skewness coefficient we use the sample 
unbiased version of the standard deviation. For the kurtosis, we do not subtract 3. 
}

\value{
If "pvalue" is FALSE, a vector with the relevant coefficient. Otherwise a matrix with two columns. 
The kurtosis or skewness coefficient and the p-value from the hypothesis test that they are significantly different
from 3 or 0 respectively.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{skew}, \link{skew.test2}, \link{colMedians}, \link{colmeans}, \link{colVars}, \link{sftests}
}
}

\examples{
## 200 variables, hence 200 F-tests will be performed
x = matrix( rnorm(200 * 50), ncol = 50 )
## 200 observations in total
colkurtosis(x)
colskewness(x)
x <- NULL
}

\keyword{ Skewness coefficient }
\keyword{ kurtosis coefficient }
