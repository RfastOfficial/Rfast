\name{MLE of continuous univariate distributions defined on the real line}
\alias{normal.mle}
\alias{gumbel.mle}
\alias{cauchy.mle}
\alias{logistic.mle}
\alias{ct.mle}
\alias{tmle}
\alias{wigner.mle}
\alias{laplace.mle}
\title{
MLE of continuous univariate distributions defined on the real line}

\description{
MLE of continuous univariate distributions defined on the real line.
}

\usage{
normal.mle(x) 
gumbel.mle(x, tol = 1e-09)
cauchy.mle(x, tol = 1e-09)
logistic.mle(x, tol = 1e-07)
ct.mle(x, tol = 1e-09)
tmle(x, v = 5, tol = 1e-08)
wigner.mle(x, tol = 1e-09)
laplace.mle(x)
}

\arguments{
\item{x}{
A numerical vector with data.
}
\item{v}{The degrees of freedom of the t distribution.}
\item{tol}{
The tolerance level up to which the maximisation stops set to 1e-09 by default.
}
}

\details{
Instead of maximising the log-likelihood via a numerical optimiser we have used a Newton-Raphson algorithm which 
is faster. See wikipedia for the equation to be solved. For the t distribution we need the degrees of freedom and 
estimate the location and scatter parameters. 

The Cauchy is the t distribution with 1 degree of freedom. If you want to fit such a distribution used the
cauchy.mle and not the t.mle with 1 degree of freedom as it's faster. The Laplace distribution is also called 
double exponential distribution.

The wigner.mle refers to the wigner semicircle distribution.
}

\value{
Usually a list with three elements, but this is not for all cases. 
\item{iters}{The number of iterations required for the Newton-Raphson to converge.}
\item{loglik}{The value of the maximised log-likelihood.}
\item{param}{The vector of the parameters.}
}

\references{
Johnson, Norman L. Kemp, Adrianne W. Kotz, Samuel (2005). Univariate Discrete Distributions (third edition). 
Hoboken, NJ: Wiley-Interscience.

https://en.wikipedia.org/wiki/Wigner_semicircle_distribution
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{zip.mle}, \link{gammamle}, \link{vm.mle}
}
}

\examples{
x <- rt(1000,10)
a <- ct.mle(x)
res<-tmle(x, v = a$nu)
res<-cauchy.mle(x)
res<-normal.mle(x)
res<-logistic.mle(x)
res<-gumbel.mle(x)
}

\keyword{ Continuous distributions }
\keyword{ t distribution }
\keyword{ Cauchy }
\keyword{ Laplace distribution }
\keyword{ Gumbel distribution }
\keyword{ Logistic distribution }
\keyword{ Wigner semicircle distribution }