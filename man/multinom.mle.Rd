\name{MLE for multivariate discrete data}
\alias{multinom.mle}
\alias{dirimultinom.mle}
\alias{colpoisson.mle}
\alias{colgeom.mle}
\title{
MLE for multivariate discrete data
}

\description{
MLE for multivariate discrete data.
}

\usage{
multinom.mle(x)
dirimultinom.mle(x, tol = 1e-07) 
colpoisson.mle(x)
colgeom.mle(x, type = 1)
}

\arguments{
\item{x}{
A matrix with discrete valued non negative data.
}
\item{tol}{the tolerance level to terminate the Newton-Raphson algorithm for the Dirichlet multinomial distribution.
}
\item{type}{
This is for the geometric distribution only. Type 1 refers to the case where the minimum is zero and type 2 for the 
case of the minimum being 1. 
}
}

\details{
For the Poisson and geometric distributions we simply fit independent Poisson and geometric distributions 
respectively.
}

\value{
A list including:
\item{loglik}{A vector with the value of the maximised log-likelihood.
}
\item{param}{A vector of the parameters.
}
}

\references{
Johnson Norman L., Kotz Samuel and Balakrishnan (1997). Discrete Multivariate Distributions. Wiley

Minka Thomas (2012). Estimating a Dirichlet distribution. Technical report. 
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{poisson.mle}, \link{zip.mle}, \link{ztp.mle}, \link{negbin.mle}, \link{poisson.nb}
}
}

\examples{
x <- t( rmultinom(1000, 20, c(0.4, 0.5, 0.1) ) )
res<-multinom.mle(x)
res<-colpoisson.mle(x)
x <- NULL
}

\keyword{ Multinomial distribution }
\keyword{ multivariate discrete data }