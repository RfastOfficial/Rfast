\name{Skewness and kurtosis coefficients}
\alias{skew}
\alias{kurt}
\title{
Skewness and kurtosis coefficients
}

\description{
Skewness and kurtosis coefficients.
}

\usage{
skew(x, pvalue = FALSE)

kurt(x, pvalue = FALSE) 
}

\arguments{
\item{x}{
A numerical vector with data.
}
\item{pvalue}{
If you want a hypothesis test that the skewness or kurtosis are significant set this to TRUE. This checks
whether the skewness is significantly different from 0 and whether the kurtosis is significantly different from
3.  
}
}

\details{
The sample skewness and kurtosis coefficient are calculated. For the kurtosis we do not subtract 3. 
} 

\value{
If "pvalue" is FALSE (default value) the skewness or kurtosis coefficients are returned. Otherwise, 
the p-value of the significance of the coefficient is returned.
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
\code{ \link{colskewness}, \link{skew.test2}, \link{colmeans}, \link{colVars}, \link{colMedians}
}
}

\examples{
x <- rgamma(500,1, 4)
res<-skew(x)
res<-kurt(x, TRUE)
}

\keyword{ Skewness }
\keyword{ kurtosis }
