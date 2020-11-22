\name{MLE of the multivariate (log-) normal distribution}
\alias{mvnorm.mle}
\alias{mvlnorm.mle}
\title{
MLE of the multivariate (log-) normal distribution
}

\description{
MLE of the multivariate (log-) normal distribution.
}

\usage{
mvnorm.mle(x)
mvlnorm.mle(x)
}

\arguments{
\item{x}{
A matrix with numerical data.
}
}

\details{
The mean vector, covariance matrix and the value of the log-likelihood of the multivariate normal or log-normal distribution is calculated. 
For the log-normal distribution we also provide the expected value and the covariance matrix. 
}

\value{
A list including:
\item{loglik}{
The log-likelihood multivariate distribution.
}
\item{mu}{
The mean vector.
}
\item{sigma}{
The covariance matrix.
}
\item{m}{
The expected mean vector of the multivariate log-normal distribution.
}
\item{s}{
The expected covariance matrix of the multivariate log-normal distribution.
}
}

\references{
Kotz, S., Balakrishnan, N., & Johnson, N. L. (2004). Continuous multivariate distributions, Volume 1: Models and applications (Vol. 1). John wiley & sons.

http://isi.cbs.nl/iamamember/CD2/pdf/329.PDF

https://en.wikipedia.org/wiki/Log-normal_distribution#Multivariate_log-normal
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{multinom.mle}, \link{dmvnorm}, \link{gaussian.nb}
}
}

\examples{
x <- matrnorm(100, 4)
res<-mvnorm.mle(x)
x <- NULL
}

\keyword{ Multivariate normal distribution }
