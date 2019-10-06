\name{MLE of count data (univariate discrete distributions)}
\alias{zip.mle}
\alias{ztp.mle}
\alias{negbin.mle}
\alias{binom.mle}
\alias{borel.mle}
\alias{geom.mle}
\alias{logseries.mle}
\alias{poisson.mle}
\alias{betageom.mle}
\alias{betabinom.mle}
\title{
MLE of count data
}

\description{
MLE of count data.
}

\usage{
zip.mle(x, tol = 1e-09)
ztp.mle(x, tol = 1e-09)
negbin.mle(x, type = 1, tol = 1e-09) 
binom.mle(x, N = NULL, tol = 1e-07)
borel.mle(x)
geom.mle(x, type = 1)
logseries.mle(x, tol = 1e-09)
poisson.mle(x)
betageom.mle(x, tol = 1e-07)
betabinom.mle(x, N, tol = 1e-07)
}

\arguments{
\item{x}{
A vector with discrete valued data.
}
\item{type}{
This argument is for the negative binomial and the geometric distribution. 
In the negative binomial you can choose which way your prefer. Type 1 is for smal sample sizes, whereas 
type 2 is for larger ones as is faster. For the geometric it is related to its two forms. Type 1 refers to the case
where the minimum is zero and type 2 for the case of the minimum being 1. 
}
\item{N}{
This is for the binomial distribution only, specifying the total number of successes. If NULL, it is sestimated by the data.
It can also be a vector of successes.
}
\item{tol}{
The tolerance level up to which the maximisation stops set to 1e-09 by default.
}
}

\details{
Instead of maximising the log-likelihood via a numerical optimiser we used a Newton-Raphson algorithm which is faster. 

See wikipedia for the equation to be solved in the case of the zero inflated distribution. https://en.wikipedia.org/wiki/Zero-inflated_model. 
In order to avoid negative values we have used link functions, log for the \eqn{lambda} and logit for the \eqn{\pi} as suggested by Lambert (1992). 
As for the zero truncated Poisson see https://en.wikipedia.org/wiki/Zero-truncated_Poisson_distribution.

zip.mle is for the zero inflated Poisson, whereas ztp.mle is for the zero truncated Poisson distribution.
}

\value{
The following list is not inclusive of all cases. Different functions have different names. In general a list including:
\item{mess}{
This is for the negbin.mle only. If there is no reason to use the negative binomial distribution a message will appear, otherwise this is NULL.
}
\item{iters}{The number of iterations required for the Newton-Raphson to converge.
}
\item{loglik}{The value of the maximised log-likelihood.
}
\item{prob}{The probability parameter of the distribution. In some distributions this argument might have a different name. 
For example, param in the zero inflated Poisson.
}
}

\references{
Lambert Diane (1992). Zero-Inflated Poisson Regression, with an Application to Defects in 
Manufacturing. Technometrics. 34 (1): 1-14

Johnson Norman L., Kotz Samuel and Kemp Adrienne W. (1992). Univariate Discrete 
Distributions (2nd ed.). Wiley
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{poisson_only}, \link{colrange}
}
}

\examples{
x <- rpois(100, 2)
res<-zip.mle(x)
res<-poisson.mle(x)
## small difference in the two log-likelihoods as expected.

x <- rpois(100, 10)
x[x == 0 ] <- 1
res<-ztp.mle(x)
res<-poisson.mle(x)
## significant difference in the two log-likelihoods. 

x <- rnbinom(100, 10, 0.6)
res<-poisson.mle(x)
res<-negbin.mle(x)
}

\keyword{ Poisson distribution }
\keyword{ zero inflated Poisson }
\keyword{ excessive zeros }
\keyword{ zero truncated Poisson }
\keyword{ negative binomial }
\keyword{ binomial distribution }
\keyword{ geometric distribution }