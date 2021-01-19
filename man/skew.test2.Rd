\name{Hypothesis testing between two skewness or kurtosis coefficients}
\alias{skew.test2}
\alias{kurt.test2}
\title{
Hypothesis testing between two skewness or kurtosis coefficients
}

\description{
Hypothesis testing between two skewness or kurtosis coefficients.
}

\usage{
skew.test2(x, y)

kurt.test2(x, y) 
}

\arguments{
\item{x}{
A numerical vector with data.
}
\item{y}{
A numerical vector with data, not necessarily of the same size.
}
}

\details{
The skewness of kurtosis coefficients between two samples are being compared.
} 

\value{
A vector with the test statistic and its associated p-value. 
}

\references{
https://en.wikipedia.org/wiki/Skewness

https://en.wikipedia.org/wiki/Kurtosis
}

\author{
Klio Lakiotaki

R implementation and documentation: Klio Lakiotaki <kliolak@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{skew}, \link{colskewness}, \link{colmeans}, \link{colVars}, \link{colMedians}
}
}

\examples{
x <- rgamma(150,1, 4)
y <- rgamma(100, 1, 4)
res<-skew.test2(x, y)
res<-kurt.test2(x, y)
}

\keyword{ Skewness }
\keyword{ kurtosis }
\keyword{ hypothesis testing }
