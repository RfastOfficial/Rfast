\name{MLE of continuous univariate distributions defined on the positive line}
\alias{gammamle}
\alias{chisq.mle}
\alias{weibull.mle}
\alias{lomax.mle}
\alias{foldnorm.mle}
\alias{betaprime.mle}
\alias{logcauchy.mle}
\alias{loglogistic.mle}
\alias{halfnorm.mle}
\alias{invgauss.mle}
\alias{lognorm.mle}
\alias{pareto.mle}
\alias{expmle}
\alias{exp2.mle}
\alias{maxboltz.mle}
\alias{rayleigh.mle}
\alias{normlog.mle}
\alias{lindley.mle}
\title{
MLE of continuous univariate distributions defined on the positive line}

\description{
MLE of continuous univariate distributions defined on the positive line.
}

\usage{
gammamle(x, tol = 1e-09) 
chisq.mle(x, tol = 1e-09)
weibull.mle(x, tol = 1e-09, maxiters = 100)
lomax.mle(x, tol = 1e-09)
foldnorm.mle(x, tol = 1e-09)
betaprime.mle(x, tol = 1e-09)
logcauchy.mle(x, tol = 1e-09)
loglogistic.mle(x, tol = 1e-09)
halfnorm.mle(x)
invgauss.mle(x)
lognorm.mle(x)
pareto.mle(x)
expmle(x)
exp2.mle(x)
maxboltz.mle(x)
rayleigh.mle(x)
normlog.mle(x)
lindley.mle(x)
}

\arguments{
\item{x}{
A vector with positive valued data (zeros are not allowed).
}
\item{tol}{
The tolerance level up to which the maximisation stops; set to 1e-09 by default.
}
\item{maxiters}{
The maximum number of iterations the Newton-Raphson will perform. 
}
}

\details{
Instead of maximising the log-likelihood via a numerical optimiser we have used a Newton-Raphson algorithm which is faster. See wikipedia for the equations to be solved. 
For the t distribution we need the degrees of freedom and estimate the location and scatter parameters. If you want to to fit an inverse gamma distribution simply do 
"gamma.mle(1/x)". The log-likelihood and the parameters are for the inverse gamma. 

The "normlog.mle" is simply the normal distribution where all values are positive. Note, this is not log-normal. It is the normal with a log link. Similarly to the inverse gaussian
distribution where the mean is an exponentiated. This comes from the GLM theory. 
 
}

\value{
Usually a list with three elements, but this is not for all cases. 
\item{iters}{The number of iterations required for the Newton-Raphson to converge.}
\item{loglik}{The value of the maximised log-likelihood.}
\item{param}{The vector of the parameters.}
}

\references{
Kalimuthu Krishnamoorthy, Meesook Lee and Wang Xiao (2015). Likelihood ratio tests for comparing several 
gamma distributions. Environmetrics, 26(8):571-583.

N.L. Johnson, S. Kotz \& N. Balakrishnan (1994). Continuous Univariate Distributions, Volume 1 (2nd Edition).

N.L. Johnson, S. Kotz \& N. Balakrishnan (1970). Distributions in statistics: continuous univariate 
distributions, Volume 2

Tsagris M., Beneki C. and Hassani H. (2014). On the folded normal distribution. Mathematics, 2(1):12-28.

Sharma V. K., Singh S. K., Singh U. \& Agiwal V. (2015). The inverse Lindley distribution: 
a stress-strength reliability model with application to head and neck cancer data. Journal of 
Industrial and Production Engineering, 32(3): 162-173.

You can also check the relevant wikipedia pages for these distributions.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{zip.mle}, \link{normal.mle}, \link{beta.mle}
}
}

\examples{
x <- rgamma(100, 3, 4)
system.time( for (i in 1:20) gammamle(x) )
## system.time( for (i in 1:20) fitdistr(x,"gamma") )
a <- glm(x ~ 1, gaussian(log) )
res<-normlog.mle(x)
}

\keyword{ Continuous distributions }
\keyword{ gamma distribution }
\keyword{ folded normal }
\keyword{ half normal }
\keyword{ Weibull }
\keyword{ Pareto }
\keyword{ beta prime }
