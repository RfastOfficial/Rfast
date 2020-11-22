\name{Multivariate kurtosis}
\alias{mvkurtosis}
\title{
Multivariate kurtosis
}

\description{
Multivariate kurtosis.
}

\usage{
mvkurtosis(x)
}

\arguments{
\item{x}{
A numerical matrix.
}
}

\details{
The multivariate kurtosis is calcualted.
} 

\value{
A number, the multivariate kurtosis.
}

\references{
K. V. Mardia (1970). Measures of Multivariate Skewness and Kurtosis with Applications
Biometrika, 57(3):519-530.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{colskewness}, \link{skew.test2}, \link{colmeans}, \link{colVars}, \link{colMedians}
}
}

\examples{
x <- as.matrix(iris[, 1:4])
res<-mvkurtosis(x)
}

\keyword{ Multivariate data }
\keyword{ kurtosis }
